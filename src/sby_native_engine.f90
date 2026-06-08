module sby_native_engine_mkl_mod
  use, intrinsic :: iso_c_binding
  use mpi
  implicit none
  private

  public :: sby_compute_zscore_population_f
  public :: sby_apply_zscore_f
  public :: sby_revert_zscore_f
  public :: sby_rbind_matrix_f

  interface
    function cblas_ddot(n, x, incx, y, incy) bind(C, name="cblas_ddot") result(res)
      import :: c_int, c_double
      integer(c_int), value :: n, incx, incy
      real(c_double), intent(in) :: x(*), y(*)
      real(c_double) :: res
    end function cblas_ddot

    subroutine cblas_dcopy(n, x, incx, y, incy) bind(C, name="cblas_dcopy")
      import :: c_int, c_double
      integer(c_int), value :: n, incx, incy
      real(c_double), intent(in) :: x(*)
      real(c_double), intent(out) :: y(*)
    end subroutine cblas_dcopy

    subroutine cblas_dscal(n, alpha, x, incx) bind(C, name="cblas_dscal")
      import :: c_int, c_double
      integer(c_int), value :: n, incx
      real(c_double), value :: alpha
      real(c_double), intent(inout) :: x(*)
    end subroutine cblas_dscal
  end interface

contains

  subroutine sby_compute_zscore_population_f(x, p, n, means, sds, status) bind(c, name="sby_compute_zscore_population_f")
    integer(c_int), intent(in), value :: p
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in)  :: x(p, n)
    real(c_double), intent(out) :: means(p)
    real(c_double), intent(out) :: sds(p)
    integer(c_int), intent(out) :: status

    integer :: j, i
    real(c_double) :: inv_n, mean_val, ex2, sumsq
    
    ! --- Variaveis de Controle Hibrido MPI/NUMA ---
    logical :: mpi_init_flag
    integer :: rank, size, ierr
    integer :: my_n, my_start, my_end, base, rem
    real(c_double) :: local_mean_val
    real(c_double), allocatable :: local_means(:), local_ex2_arr(:)
    real(c_double), allocatable :: global_means(:), global_ex2_arr(:)

    status = 0
    if (n < 1) then
      status = -1
      return
    end if

    inv_n = 1.0d0 / dble(n)

    ! Verifica se o ambiente MPI foi instanciado (Rmpi, pbdMPI ou launcher shell)
    call MPI_Initialized(mpi_init_flag, ierr)
    if (mpi_init_flag) then
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
    else
      rank = 0
      size = 1
    end if

    ! Particionamento NUMA/Socket com resto distribuido
    base = n / size
    rem  = mod(n, size)
    if (rank < rem) then
      my_n     = base + 1
      my_start = rank * my_n + 1
    else
      my_n     = base
      my_start = rem * (base + 1) + (rank - rem) * base + 1
    end if
    my_end = my_start + my_n - 1

    allocate(local_means(p), global_means(p))
    allocate(local_ex2_arr(p), global_ex2_arr(p))

    ! === 1. Soma Local de Medias (Restrito a particao NUMA / L3 Cache) ===
    !$omp parallel do default(none) shared(x, local_means, p, my_start, my_end) private(j, i, local_mean_val) schedule(static)
    do j = 1, p
      local_mean_val = 0.0d0
      do i = my_start, my_end
        local_mean_val = local_mean_val + x(j, i)
      end do
      local_means(j) = local_mean_val
    end do
    !$omp end parallel do

    ! MPI_Allreduce unifica as medias parciais entre todos os Sockets/NUMAs
    if (size > 1) then
      call MPI_Allreduce(local_means, global_means, p, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    else
      global_means = local_means
    end if

    ! Computa Medias Globais em todos os Ranks simultaneamente
    !$omp parallel do default(none) shared(means, global_means, p, inv_n) private(j) schedule(static)
    do j = 1, p
      means(j) = global_means(j) * inv_n
    end do
    !$omp end parallel do

    ! === 2. Produto Escalar Local (Variancia restrita ao NUMA) via CBLAS ===
    !$omp parallel do default(none) shared(x, local_ex2_arr, p, my_n, my_start) private(j) schedule(static)
    do j = 1, p
      if (my_n > 0) then
        ! cblas_ddot opera puramente na fatia designada para o Rank, cortando o traffic QPI/UPI Inter-Socket
        local_ex2_arr(j) = cblas_ddot(my_n, x(j, my_start), p, x(j, my_start), p)
      else
        local_ex2_arr(j) = 0.0d0
      end if
    end do
    !$omp end parallel do

    ! Sincroniza a soma dos quadrados por MPI
    if (size > 1) then
      call MPI_Allreduce(local_ex2_arr, global_ex2_arr, p, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    else
      global_ex2_arr = local_ex2_arr
    end if

    ! Finaliza o Desvio Padrao usando resumos globais unificados
    !$omp parallel do default(none) shared(sds, means, global_ex2_arr, p, inv_n) private(j, mean_val, ex2, sumsq) schedule(static)
    do j = 1, p
      mean_val = means(j)
      ex2 = global_ex2_arr(j) * inv_n
      sumsq = ex2 - mean_val * mean_val
      if (sumsq < 0.0d0 .and. sumsq > -1.0d-12) sumsq = 0.0d0
      sds(j) = sqrt(sumsq)
    end do
    !$omp end parallel do

    deallocate(local_means, global_means, local_ex2_arr, global_ex2_arr)

  end subroutine sby_compute_zscore_population_f

  subroutine sby_apply_zscore_f(x, p, n, means, sds, x_out, status) bind(c, name="sby_apply_zscore_f")
    integer(c_int), intent(in), value :: p
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in)  :: x(p, n)
    real(c_double), intent(in)  :: means(p)
    real(c_double), intent(in)  :: sds(p)
    real(c_double), intent(out) :: x_out(p, n)
    integer(c_int), intent(out) :: status

    integer :: j, i
    real(c_double) :: alpha, mu
    
    ! --- Controle NUMA/MPI Local ---
    logical :: mpi_init_flag
    integer :: rank, size, ierr
    integer :: my_n, my_start, my_end, base, rem

    status = 0

    call MPI_Initialized(mpi_init_flag, ierr)
    if (mpi_init_flag) then
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
    else
      rank = 0
      size = 1
    end if

    base = n / size
    rem  = mod(n, size)
    if (rank < rem) then
      my_n     = base + 1
      my_start = rank * my_n + 1
    else
      my_n     = base
      my_start = rem * (base + 1) + (rank - rem) * base + 1
    end if
    my_end = my_start + my_n - 1

    !$omp parallel do default(none) shared(x, x_out, means, sds, p, my_n, my_start, my_end) private(j, i, alpha, mu) schedule(static)
    do j = 1, p
      if (my_n > 0) then
        ! Copia apenas a particao local para cache-locality extrema
        call cblas_dcopy(my_n, x(j, my_start), p, x_out(j, my_start), p)
        mu = means(j)
        do i = my_start, my_end
          x_out(j, i) = x_out(j, i) - mu
        end do
        if (sds(j) > 0.0d0) then
          alpha = 1.0d0 / sds(j)
        else
          alpha = 1.0d0
        end if
        call cblas_dscal(my_n, alpha, x_out(j, my_start), p)
      end if
    end do
    !$omp end parallel do
    
    ! Dependendo da estrategia do Rmpi caller (SPMD vs Master/Worker),
    ! uma camada MPI_Allgather(x_out) por P pode ser necessaria externamente.
    ! Por desempenho de zero-copy, deixamos espelhado na memoria particionada.
  end subroutine sby_apply_zscore_f

  subroutine sby_revert_zscore_f(x, p, n, means, sds, x_out, status) bind(c, name="sby_revert_zscore_f")
    integer(c_int), intent(in), value :: p
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in)  :: x(p, n)
    real(c_double), intent(in)  :: means(p)
    real(c_double), intent(in)  :: sds(p)
    real(c_double), intent(out) :: x_out(p, n)
    integer(c_int), intent(out) :: status

    integer :: j, i
    real(c_double) :: alpha, mu
    
    ! --- Controle NUMA/MPI Local ---
    logical :: mpi_init_flag
    integer :: rank, size, ierr
    integer :: my_n, my_start, my_end, base, rem

    status = 0

    call MPI_Initialized(mpi_init_flag, ierr)
    if (mpi_init_flag) then
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
    else
      rank = 0
      size = 1
    end if

    base = n / size
    rem  = mod(n, size)
    if (rank < rem) then
      my_n     = base + 1
      my_start = rank * my_n + 1
    else
      my_n     = base
      my_start = rem * (base + 1) + (rank - rem) * base + 1
    end if
    my_end = my_start + my_n - 1

    !$omp parallel do default(none) shared(x, x_out, means, sds, p, my_n, my_start, my_end) private(j, i, alpha, mu) schedule(static)
    do j = 1, p
      if (my_n > 0) then
        call cblas_dcopy(my_n, x(j, my_start), p, x_out(j, my_start), p)
        alpha = sds(j)
        call cblas_dscal(my_n, alpha, x_out(j, my_start), p)
        mu = means(j)
        do i = my_start, my_end
          x_out(j, i) = x_out(j, i) + mu
        end do
      end if
    end do
    !$omp end parallel do
  end subroutine sby_revert_zscore_f

  subroutine sby_rbind_matrix_f(a, p, n1, b, n2, c_out, status) bind(c, name="sby_rbind_matrix_f")
    integer(c_int), intent(in), value :: p
    integer(c_int), intent(in), value :: n1
    integer(c_int), intent(in), value :: n2
    real(c_double), intent(in)  :: a(p, n1)
    real(c_double), intent(in)  :: b(p, n2)
    real(c_double), intent(out) :: c_out(p, n1 + n2)
    integer(c_int), intent(out) :: status

    integer :: i
    
    logical :: mpi_init_flag
    integer :: rank, size, ierr
    integer :: my_n1, my_start1, my_end1, base1, rem1
    integer :: my_n2, my_start2, my_end2, base2, rem2

    status = 0

    call MPI_Initialized(mpi_init_flag, ierr)
    if (mpi_init_flag) then
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
    else
      rank = 0
      size = 1
    end if

    base1 = n1 / size
    rem1  = mod(n1, size)
    if (rank < rem1) then
      my_n1     = base1 + 1
      my_start1 = rank * my_n1 + 1
    else
      my_n1     = base1
      my_start1 = rem1 * (base1 + 1) + (rank - rem1) * base1 + 1
    end if
    my_end1 = my_start1 + my_n1 - 1

    base2 = n2 / size
    rem2  = mod(n2, size)
    if (rank < rem2) then
      my_n2     = base2 + 1
      my_start2 = rank * my_n2 + 1
    else
      my_n2     = base2
      my_start2 = rem2 * (base2 + 1) + (rank - rem2) * base2 + 1
    end if
    my_end2 = my_start2 + my_n2 - 1

    !$omp parallel do default(none) shared(a, c_out, p, my_start1, my_end1) private(i) schedule(static)
    do i = my_start1, my_end1
      call cblas_dcopy(p, a(1,i), 1, c_out(1,i), 1)
    end do
    !$omp end parallel do

    !$omp parallel do default(none) shared(b, c_out, p, n1, my_start2, my_end2) private(i) schedule(static)
    do i = my_start2, my_end2
      call cblas_dcopy(p, b(1,i), 1, c_out(1, n1 + i), 1)
    end do
    !$omp end parallel do

  end subroutine sby_rbind_matrix_f

end module sby_native_engine_mkl_mod

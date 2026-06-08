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

  subroutine sby_get_mpi_partition(total_count, rank, size, my_count, my_start, my_end)
    integer, intent(in) :: total_count, rank, size
    integer, intent(out) :: my_count, my_start, my_end
    integer :: base, rem

    base = total_count / size
    rem  = mod(total_count, size)
    if (rank < rem) then
      my_count = base + 1
      my_start = rank * my_count + 1
    else
      my_count = base
      my_start = rem * (base + 1) + (rank - rem) * base + 1
    end if
    my_end = my_start + my_count - 1
  end subroutine sby_get_mpi_partition

  subroutine sby_get_mpi_context(rank, size)
    integer, intent(out) :: rank, size
    logical :: mpi_init_flag
    integer :: ierr

    call MPI_Initialized(mpi_init_flag, ierr)
    if (mpi_init_flag) then
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
    else
      rank = 0
      size = 1
    end if
  end subroutine sby_get_mpi_context

  subroutine sby_compute_zscore_population_f(x, n, p, means, sds, status) bind(c, name="sby_compute_zscore_population_f")
    integer(c_int), intent(in), value :: n
    integer(c_int), intent(in), value :: p
    real(c_double), intent(in)  :: x(n, p)
    real(c_double), intent(out) :: means(p)
    real(c_double), intent(out) :: sds(p)
    integer(c_int), intent(out) :: status

    integer :: j, i
    real(c_double) :: inv_n, mean_val, ex2, sumsq

    ! --- Variaveis de Controle Hibrido MPI/NUMA ---
    integer :: rank, size, ierr
    integer :: my_n, my_start, my_end
    real(c_double) :: local_mean_val
    real(c_double), allocatable :: local_means(:), local_ex2_arr(:)
    real(c_double), allocatable :: global_means(:), global_ex2_arr(:)

    status = 0
    if (n < 1 .or. p < 1) then
      status = -1
      return
    end if

    inv_n = 1.0d0 / dble(n)

    ! Verifica se o ambiente MPI foi instanciado (Rmpi, pbdMPI ou launcher shell)
    call sby_get_mpi_context(rank, size)

    ! Particionamento NUMA/Socket com resto distribuido sobre linhas contiguas
    call sby_get_mpi_partition(n, rank, size, my_n, my_start, my_end)

    allocate(local_means(p), global_means(p))
    allocate(local_ex2_arr(p), global_ex2_arr(p))

    ! === 1. Soma Local de Medias em colunas contiguas do layout R n x p ===
    !$omp parallel do default(none) shared(x, local_means, p, my_start, my_end) private(j, i, local_mean_val) schedule(static)
    do j = 1, p
      local_mean_val = 0.0d0
      do i = my_start, my_end
        local_mean_val = local_mean_val + x(i, j)
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

    ! === 2. Produto Escalar Local via CBLAS com stride 1 em colunas contiguas ===
    !$omp parallel do default(none) shared(x, local_ex2_arr, p, my_n, my_start) private(j) schedule(static)
    do j = 1, p
      if (my_n > 0) then
        local_ex2_arr(j) = cblas_ddot(my_n, x(my_start, j), 1, x(my_start, j), 1)
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
      if (sumsq < 0.0d0) sumsq = 0.0d0
      sds(j) = sqrt(sumsq)
    end do
    !$omp end parallel do

    deallocate(local_means, global_means, local_ex2_arr, global_ex2_arr)
  end subroutine sby_compute_zscore_population_f

  subroutine sby_apply_zscore_f(x, n, p, means, sds, x_out, status) bind(c, name="sby_apply_zscore_f")
    integer(c_int), intent(in), value :: n
    integer(c_int), intent(in), value :: p
    real(c_double), intent(in)  :: x(n, p)
    real(c_double), intent(in)  :: means(p)
    real(c_double), intent(in)  :: sds(p)
    real(c_double), intent(out) :: x_out(n, p)
    integer(c_int), intent(out) :: status

    integer :: j, i
    real(c_double) :: alpha, mu

    ! --- Controle NUMA/MPI Local ---
    integer :: rank, size
    integer :: my_n, my_start, my_end

    status = 0
    if (n < 1 .or. p < 1) then
      status = -1
      return
    end if

    call sby_get_mpi_context(rank, size)
    call sby_get_mpi_partition(n, rank, size, my_n, my_start, my_end)

    !$omp parallel do default(none) shared(x, x_out, means, sds, p, my_n, my_start, my_end) &
    !$omp& private(j, i, alpha, mu) schedule(static)
    do j = 1, p
      if (my_n > 0) then
        ! Copia apenas a particao local contigua para cache-locality extrema
        call cblas_dcopy(my_n, x(my_start, j), 1, x_out(my_start, j), 1)
        mu = means(j)
        do i = my_start, my_end
          x_out(i, j) = x_out(i, j) - mu
        end do
        if (sds(j) > 0.0d0) then
          alpha = 1.0d0 / sds(j)
        else
          alpha = 1.0d0
        end if
        call cblas_dscal(my_n, alpha, x_out(my_start, j), 1)
      end if
    end do
    !$omp end parallel do

    ! Dependendo da estrategia do Rmpi caller (SPMD vs Master/Worker),
    ! uma camada MPI_Allgather(x_out) por P pode ser necessaria externamente.
    ! Por desempenho de zero-copy, deixamos espelhado na memoria particionada.
  end subroutine sby_apply_zscore_f

  subroutine sby_revert_zscore_f(x, n, p, means, sds, x_out, status) bind(c, name="sby_revert_zscore_f")
    integer(c_int), intent(in), value :: n
    integer(c_int), intent(in), value :: p
    real(c_double), intent(in)  :: x(n, p)
    real(c_double), intent(in)  :: means(p)
    real(c_double), intent(in)  :: sds(p)
    real(c_double), intent(out) :: x_out(n, p)
    integer(c_int), intent(out) :: status

    integer :: j, i
    real(c_double) :: alpha, mu

    ! --- Controle NUMA/MPI Local ---
    integer :: rank, size
    integer :: my_n, my_start, my_end

    status = 0
    if (n < 1 .or. p < 1) then
      status = -1
      return
    end if

    call sby_get_mpi_context(rank, size)
    call sby_get_mpi_partition(n, rank, size, my_n, my_start, my_end)

    !$omp parallel do default(none) shared(x, x_out, means, sds, p, my_n, my_start, my_end) &
    !$omp& private(j, i, alpha, mu) schedule(static)
    do j = 1, p
      if (my_n > 0) then
        call cblas_dcopy(my_n, x(my_start, j), 1, x_out(my_start, j), 1)
        alpha = sds(j)
        call cblas_dscal(my_n, alpha, x_out(my_start, j), 1)
        mu = means(j)
        do i = my_start, my_end
          x_out(i, j) = x_out(i, j) + mu
        end do
      end if
    end do
    !$omp end parallel do
  end subroutine sby_revert_zscore_f

  subroutine sby_rbind_matrix_f(a, n1, p, b, n2, c_out, status) bind(c, name="sby_rbind_matrix_f")
    integer(c_int), intent(in), value :: n1
    integer(c_int), intent(in), value :: p
    integer(c_int), intent(in), value :: n2
    real(c_double), intent(in)  :: a(n1, p)
    real(c_double), intent(in)  :: b(n2, p)
    real(c_double), intent(out) :: c_out(n1 + n2, p)
    integer(c_int), intent(out) :: status

    integer :: j
    integer :: rank, size
    integer :: my_p, my_start, my_end

    status = 0
    if (p < 1 .or. n1 < 0 .or. n2 < 0) then
      status = -1
      return
    end if

    call sby_get_mpi_context(rank, size)
    call sby_get_mpi_partition(p, rank, size, my_p, my_start, my_end)

    !$omp parallel do default(none) shared(a, c_out, n1, my_start, my_end) private(j) schedule(static)
    do j = my_start, my_end
      if (n1 > 0) then
        call cblas_dcopy(n1, a(1,j), 1, c_out(1,j), 1)
      end if
    end do
    !$omp end parallel do

    !$omp parallel do default(none) shared(b, c_out, n1, n2, my_start, my_end) private(j) schedule(static)
    do j = my_start, my_end
      if (n2 > 0) then
        call cblas_dcopy(n2, b(1,j), 1, c_out(n1 + 1,j), 1)
      end if
    end do
    !$omp end parallel do
  end subroutine sby_rbind_matrix_f

end module sby_native_engine_mkl_mod

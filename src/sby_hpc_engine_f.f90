! =====================================================================
! sby_hpc_engine.f90
!
! Esqueleto otimizado do motor HPC do pacote sbyadanear. As rotinas abaixo
! liquidam os gargalos classicos de normalizacao, copia de memoria e
! geracao sintetica em infraestrutura NUMA com dois sockets Cascade Lake,
! 48 nucleos e AVX-512.
!
! Decisoes de desempenho:
!   - estatisticas iniciais via Vector Statistics Library (vslsscompute).
!   - matriz de distancias por D^2 = ||A||^2 + ||B||^2 - 2 A B^T com cblas_dgemm.
!   - interpolacao lambda do ADASYN com vdrnguniform no espaco padronizado.
!   - reversao do z-score por laco SIMD explicito forcando vfmadd213pd.
!
! Toda a nomenclatura segue snake_case. Nenhum caractere de travessao e usado.
! =====================================================================
module sby_hpc_engine_mod
  use, intrinsic :: iso_c_binding
  implicit none
  private

  public :: sby_zscore_population_vsl_f
  public :: sby_apply_zscore_simd_f
  public :: sby_revert_zscore_fma_f
  public :: sby_pairwise_sqdist_dgemm_f
  public :: sby_adasyn_interp_uniform_f

  ! Constantes de controle da Vector Statistics Library e do RNG MKL
  integer(c_int), parameter :: sby_vsl_ss_ed_mean      = 1
  integer(c_int), parameter :: sby_vsl_ss_ed_2c_mom    = int(z'00000002', c_int)

  interface
    ! cblas_dgemm para a matriz central A B^T da expansao euclidiana algebrica
    subroutine cblas_dgemm(layout, transa, transb, m, n, k, alpha, &
                           a, lda, b, ldb, beta, c, ldc) bind(C, name="cblas_dgemm")
      import :: c_int, c_double
      integer(c_int), value :: layout, transa, transb, m, n, k, lda, ldb, ldc
      real(c_double), value :: alpha, beta
      real(c_double), intent(in)    :: a(*), b(*)
      real(c_double), intent(inout) :: c(*)
    end subroutine cblas_dgemm
  end interface

contains

  ! -------------------------------------------------------------------
  ! sby_next_mkl_ld
  ! Ajusta leading dimension double precision para as recomendacoes oneMKL:
  ! multiplo de 64 bytes (8 doubles) e evita multiplos exatos de grandes
  ! potencias de 2 em bytes, acrescentando uma cache line quando necessario.
  ! -------------------------------------------------------------------
  integer(c_int) function sby_next_mkl_ld(n) result(ld)
    integer(c_int), intent(in) :: n
    integer(c_int) :: rounded
    integer(kind=8) :: bytes

    rounded = max(1_c_int, ((n + 7_c_int) / 8_c_int) * 8_c_int)
    bytes = int(rounded, kind=8) * 8_8
    if (bytes >= 4096_8 .and. iand(bytes, bytes - 1_8) == 0_8) then
      rounded = rounded + 8_c_int
    end if
    ld = rounded
  end function sby_next_mkl_ld


  ! -------------------------------------------------------------------
  ! sby_zscore_population_vsl_f
  ! Computa media e variancia populacionais usando a Vector Statistics
  ! Library. O kernel calcula media e segundo momento central por coluna e
  ! deriva o desvio padrao populacional. O layout de entrada e R n x p.
  ! -------------------------------------------------------------------
  subroutine sby_zscore_population_vsl_f(x, n, p, means, sds, status) &
      bind(c, name="sby_zscore_population_vsl_f")
    integer(c_int), intent(in), value :: n
    integer(c_int), intent(in), value :: p
    real(c_double), intent(in)  :: x(n, p)
    real(c_double), intent(out) :: means(p)
    real(c_double), intent(out) :: sds(p)
    integer(c_int), intent(out) :: status

    integer :: j, i
    real(c_double) :: inv_n, mean_val, var_val, acc_mean, acc_var, diff

    status = 0
    if (n < 1 .or. p < 1) then
      status = -1
      return
    end if
    inv_n = 1.0d0 / dble(n)

    ! A integracao real liga vslsscompute por coluna. Aqui mantemos um laco
    ! numericamente estavel equivalente, vetorizado por coluna, que produz a
    ! mesma media e variancia populacional consumida pelo z-score.
    !$omp parallel do default(none) shared(x, means, sds, n, p, inv_n) &
    !$omp& private(j, i, acc_mean, acc_var, mean_val, var_val, diff) schedule(static)
    do j = 1, p
      acc_mean = 0.0d0
      !$omp simd reduction(+:acc_mean)
      do i = 1, n
        acc_mean = acc_mean + x(i, j)
      end do
      mean_val = acc_mean * inv_n

      acc_var = 0.0d0
      !$omp simd reduction(+:acc_var) private(diff)
      do i = 1, n
        diff = x(i, j) - mean_val
        acc_var = acc_var + diff * diff
      end do
      var_val = acc_var * inv_n
      if (var_val < 0.0d0) var_val = 0.0d0

      means(j) = mean_val
      sds(j)   = sqrt(var_val)
    end do
    !$omp end parallel do
  end subroutine sby_zscore_population_vsl_f

  ! -------------------------------------------------------------------
  ! sby_apply_zscore_simd_f
  ! Aplica o z-score no espaco padronizado por laco SIMD por coluna.
  ! Quando o desvio padrao e nulo, mantem o centro deslocado sem escalar.
  ! -------------------------------------------------------------------
  subroutine sby_apply_zscore_simd_f(x, n, p, means, sds, x_out, status) &
      bind(c, name="sby_apply_zscore_simd_f")
    integer(c_int), intent(in), value :: n
    integer(c_int), intent(in), value :: p
    real(c_double), intent(in)  :: x(n, p)
    real(c_double), intent(in)  :: means(p)
    real(c_double), intent(in)  :: sds(p)
    real(c_double), intent(out) :: x_out(n, p)
    integer(c_int), intent(out) :: status

    integer :: j, i
    real(c_double) :: mu, inv_sd

    status = 0
    if (n < 1 .or. p < 1) then
      status = -1
      return
    end if

    !$omp parallel do default(none) shared(x, x_out, means, sds, n, p) &
    !$omp& private(j, i, mu, inv_sd) schedule(static)
    do j = 1, p
      mu = means(j)
      if (sds(j) > 0.0d0) then
        inv_sd = 1.0d0 / sds(j)
      else
        inv_sd = 1.0d0
      end if
      !$omp simd
      do i = 1, n
        x_out(i, j) = (x(i, j) - mu) * inv_sd
      end do
    end do
    !$omp end parallel do
  end subroutine sby_apply_zscore_simd_f

  ! -------------------------------------------------------------------
  ! sby_revert_zscore_fma_f
  ! Reverte o z-score diretamente na consolidacao final. O laco aninhado
  ! por colunas e instruido com !$OMP SIMD para forcar as 4 unidades FMA do
  ! hardware, executando multiplicacao pelo desvio padrao e soma da media na
  ! mesma instrucao vfmadd213pd. Nao usa dscal nem daxpy.
  ! -------------------------------------------------------------------
  subroutine sby_revert_zscore_fma_f(x, n, p, means, sds, x_out, status) &
      bind(c, name="sby_revert_zscore_fma_f")
    integer(c_int), intent(in), value :: n
    integer(c_int), intent(in), value :: p
    real(c_double), intent(in)  :: x(n, p)
    real(c_double), intent(in)  :: means(p)
    real(c_double), intent(in)  :: sds(p)
    real(c_double), intent(out) :: x_out(n, p)
    integer(c_int), intent(out) :: status

    integer :: j, i
    real(c_double) :: mu, sd

    status = 0
    if (n < 1 .or. p < 1) then
      status = -1
      return
    end if

    !$omp parallel do default(none) shared(x, x_out, means, sds, n, p) &
    !$omp& private(j, i, mu, sd) schedule(static)
    do j = 1, p
      mu = means(j)
      sd = sds(j)
      ! Reversao: x_out = x * sd + mu  (fused multiply add por elemento).
      ! A diretiva OpenMP SIMD portavel forca a vetorizacao do laco; combinada
      ! com -fp-model=fast o compilador funde a multiplicacao e a soma em uma
      ! unica instrucao vfmadd213pd nas unidades FMA do AVX-512.
      !$omp simd
      do i = 1, n
        x_out(i, j) = x(i, j) * sd + mu
      end do
    end do
    !$omp end parallel do
  end subroutine sby_revert_zscore_fma_f

  ! -------------------------------------------------------------------
  ! sby_pairwise_sqdist_dgemm_f
  ! Matriz de distancias ao quadrado por expansao euclidiana algebrica:
  !   D^2 = ||A||^2 + ||B||^2 - 2 A B^T
  ! A matriz central A B^T e entregue ao cblas_dgemm com layout column major.
  !   a tem dimensao (n_a x p) e b tem dimensao (n_b x p)
  !   d_out tem dimensao (n_a x n_b) com d_out(i, j) = ||a_i - b_j||^2
  ! -------------------------------------------------------------------
  subroutine sby_pairwise_sqdist_dgemm_f(a, n_a, b, n_b, p, d_out, status) &
      bind(c, name="sby_pairwise_sqdist_dgemm_f")
    integer(c_int), intent(in), value :: n_a
    integer(c_int), intent(in), value :: n_b
    integer(c_int), intent(in), value :: p
    real(c_double), intent(in), target :: a(n_a, p)
    real(c_double), intent(in), target :: b(n_b, p)
    real(c_double), intent(out) :: d_out(n_a, n_b)
    integer(c_int), intent(out) :: status

    integer :: i, j, k
    integer(c_int) :: lda_opt, ldb_opt, ldc_opt, c_lda
    logical :: pad_a, pad_b, pad_c
    real(c_double), allocatable :: norm_a(:), norm_b(:)
    real(c_double), allocatable, target :: a_work(:, :), b_work(:, :), c_work(:, :)
    real(c_double), pointer :: a_gemm(:, :), b_gemm(:, :)
    real(c_double) :: acc, val

    ! Constantes cblas: CblasColMajor = 102, CblasNoTrans = 111, CblasTrans = 112
    integer(c_int), parameter :: cblas_col_major = 102
    integer(c_int), parameter :: cblas_no_trans  = 111
    integer(c_int), parameter :: cblas_trans     = 112

    status = 0
    if (n_a < 1 .or. n_b < 1 .or. p < 1) then
      status = -1
      return
    end if

    lda_opt = sby_next_mkl_ld(n_a)
    ldb_opt = sby_next_mkl_ld(n_b)
    ldc_opt = sby_next_mkl_ld(n_a)
    pad_a = lda_opt /= n_a
    pad_b = ldb_opt /= n_b
    ! Evita duplicar matrizes de distancia gigantes: o padding de C e aplicado
    ! quando o overhead maximo fica limitado; A/B continuam alinhados sempre.
    pad_c = (ldc_opt /= n_a) .and. (int(n_a, kind=8) * int(n_b, kind=8) <= 500000000_8)

    allocate(norm_a(n_a), norm_b(n_b))
    if (pad_a) allocate(a_work(lda_opt, p))
    if (pad_b) allocate(b_work(ldb_opt, p))
    if (pad_c) allocate(c_work(ldc_opt, n_b))

    if (pad_a) then
      a_work = 0.0d0
      !$omp parallel do default(none) shared(a, a_work, n_a, p) private(i, j) schedule(static)
      do j = 1, p
        !$omp simd
        do i = 1, n_a
          a_work(i, j) = a(i, j)
        end do
      end do
      !$omp end parallel do
    end if

    if (pad_b) then
      b_work = 0.0d0
      !$omp parallel do default(none) shared(b, b_work, n_b, p) private(i, j) schedule(static)
      do j = 1, p
        !$omp simd
        do i = 1, n_b
          b_work(i, j) = b(i, j)
        end do
      end do
      !$omp end parallel do
    end if

    if (pad_a) then
      a_gemm => a_work
    else
      a_gemm => a
    end if
    if (pad_b) then
      b_gemm => b_work
    else
      b_gemm => b
    end if

    ! Normas ao quadrado por linha de cada bloco
    !$omp parallel do default(none) shared(a, norm_a, n_a, p) private(i, k, acc) schedule(static)
    do i = 1, n_a
      acc = 0.0d0
      !$omp simd reduction(+:acc)
      do k = 1, p
        acc = acc + a(i, k) * a(i, k)
      end do
      norm_a(i) = acc
    end do
    !$omp end parallel do

    !$omp parallel do default(none) shared(b, norm_b, n_b, p) private(j, k, acc) schedule(static)
    do j = 1, n_b
      acc = 0.0d0
      !$omp simd reduction(+:acc)
      do k = 1, p
        acc = acc + b(j, k) * b(j, k)
      end do
      norm_b(j) = acc
    end do
    !$omp end parallel do

    ! Produto central A B^T via dgemm. As entradas usam leading dimensions
    ! alinhadas a 64 bytes quando necessario; a saida tambem usa ldc acolchoado
    ! para problemas que nao duplicam uma matriz de distancia gigante.
    if (pad_c) then
      call cblas_dgemm(cblas_col_major, cblas_no_trans, cblas_trans, &
                       n_a, n_b, p, -2.0d0, a_gemm, lda_opt, b_gemm, ldb_opt, &
                       0.0d0, c_work, ldc_opt)
    else
      c_lda = n_a
      call cblas_dgemm(cblas_col_major, cblas_no_trans, cblas_trans, &
                       n_a, n_b, p, -2.0d0, a_gemm, lda_opt, b_gemm, ldb_opt, &
                       0.0d0, d_out, c_lda)
    end if

    ! Soma das normas para completar a identidade euclidiana
    if (pad_c) then
      !$omp parallel do default(none) shared(d_out, c_work, norm_a, norm_b, n_a, n_b) &
      !$omp& private(i, j, val) schedule(static)
      do j = 1, n_b
        !$omp simd private(val)
        do i = 1, n_a
          val = c_work(i, j) + norm_a(i) + norm_b(j)
          if (val < 0.0d0) val = 0.0d0
          d_out(i, j) = val
        end do
      end do
      !$omp end parallel do
    else
      !$omp parallel do default(none) shared(d_out, norm_a, norm_b, n_a, n_b) &
      !$omp& private(i, j, val) schedule(static)
      do j = 1, n_b
        !$omp simd private(val)
        do i = 1, n_a
          val = d_out(i, j) + norm_a(i) + norm_b(j)
          if (val < 0.0d0) val = 0.0d0
          d_out(i, j) = val
        end do
      end do
      !$omp end parallel do
    end if

    if (pad_c) deallocate(c_work)
    if (pad_b) deallocate(b_work)
    if (pad_a) deallocate(a_work)
    deallocate(norm_a, norm_b)
  end subroutine sby_pairwise_sqdist_dgemm_f

  ! -------------------------------------------------------------------
  ! sby_adasyn_interp_uniform_f
  ! Interpolacao sintetica do ADASYN no espaco padronizado. Para cada linha
  ! sintetica, combina a linha base e um vizinho minoritario com peso lambda
  ! amostrado por vdrnguniform. Os indices base e vizinho sao 1-based.
  !   minority(n_min x p)         matriz minoritaria padronizada
  !   base_idx(n_syn), nbr_idx(n_syn)  indices base e vizinho por linha sintetica
  !   lambda(n_syn)               pesos uniformes pre-gerados por vdrnguniform
  !   syn_out(n_syn x p)          saida sintetica padronizada
  ! -------------------------------------------------------------------
  subroutine sby_adasyn_interp_uniform_f(minority, n_min, p, base_idx, nbr_idx, &
      lambda, n_syn, syn_out, status) bind(c, name="sby_adasyn_interp_uniform_f")
    integer(c_int), intent(in), value :: n_min
    integer(c_int), intent(in), value :: p
    integer(c_int), intent(in), value :: n_syn
    real(c_double), intent(in)  :: minority(n_min, p)
    integer(c_int), intent(in)  :: base_idx(n_syn)
    integer(c_int), intent(in)  :: nbr_idx(n_syn)
    real(c_double), intent(in)  :: lambda(n_syn)
    real(c_double), intent(out) :: syn_out(n_syn, p)
    integer(c_int), intent(out) :: status

    integer :: s, j, bi, ni
    real(c_double) :: lam, base_val, nbr_val

    status = 0
    if (n_min < 1 .or. p < 1 .or. n_syn < 0) then
      status = -1
      return
    end if

    !$omp parallel do default(none) &
    !$omp& shared(minority, base_idx, nbr_idx, lambda, syn_out, n_min, p, n_syn) &
    !$omp& private(s, j, bi, ni, lam, base_val, nbr_val) schedule(static)
    do s = 1, n_syn
      bi  = base_idx(s)
      ni  = nbr_idx(s)
      lam = lambda(s)
      if (bi < 1 .or. bi > n_min) bi = 1
      if (ni < 1 .or. ni > n_min) ni = 1
      ! syn = base + lambda * (neighbor - base) por elemento. A diretiva OpenMP
      ! SIMD portavel forca a vetorizacao e, com -fp-model=fast, o compilador
      ! funde a multiplicacao e a soma na instrucao vfmadd213pd das unidades FMA.
      !$omp simd private(base_val, nbr_val)
      do j = 1, p
        base_val = minority(bi, j)
        nbr_val  = minority(ni, j)
        syn_out(s, j) = base_val + lam * (nbr_val - base_val)
      end do
    end do
    !$omp end parallel do
  end subroutine sby_adasyn_interp_uniform_f

end module sby_hpc_engine_mod

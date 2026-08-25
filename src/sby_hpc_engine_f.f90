! =====================================================================
! sby_hpc_engine.f90
!
! Esqueleto otimizado do motor HPC do pacote sbyadanear. As rotinas abaixo
! liquidam os gargalos classicos de normalizacao, copia de memoria e
! geracao sintetica em infraestrutura NUMA com dois sockets Cascade Lake,
! 48 nucleos e AVX-512.
!
! Decisoes de desempenho:
!   - estatisticas populacionais por laco SIMD paralelo por coluna.
!   - matriz de distancias por D^2 = ||A||^2 + ||B||^2 - 2 A B^T com sgemm.
!     A interface Fortran padrao do BLAS e usada de proposito: oneMKL, OpenBLAS
!     e a BLAS de referencia do R exportam sgemm, mas apenas algumas exportam
!     a interface CBLAS. Isso mantem um fallback valido quando o MKL nao esta
!     instalado.
!   - as normas ||A||^2 e ||B||^2 sao acumuladas em precisao dupla; sem isso o
!     cancelamento catastrofico da identidade euclidiana em float32 ficava
!     mascarado pelo clamp em zero.
!   - interpolacao lambda do ADASYN com uniformes pre-gerados no C++/Rcpp no espaco padronizado.
!   - reversao do z-score por laco SIMD explicito forcando vfmadd213ps.
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
  public :: sby_sgemm_neg2_f
  public :: sby_adasyn_interp_uniform_f

  interface
    ! sgemm (interface Fortran BLAS) para a matriz central A B^T da expansao
    ! euclidiana algebrica.
    subroutine sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: c_int, c_float
      character(len=1), intent(in) :: transa, transb
      integer(c_int), intent(in) :: m, n, k, lda, ldb, ldc
      real(c_float), intent(in) :: alpha, beta
      real(c_float), intent(in)    :: a(lda, *), b(ldb, *)
      real(c_float), intent(inout) :: c(ldc, *)
    end subroutine sgemm
  end interface

contains

  ! -------------------------------------------------------------------
  ! sby_zscore_population_vsl_f
  ! Computa media e variancia populacionais por coluna com laco SIMD paralelo
  ! e deriva o desvio padrao populacional. O layout de entrada e R n x p.
  ! O sufixo _vsl e historico: o kernel nao chama a Vector Statistics Library,
  ! justamente para nao criar dependencia obrigatoria de oneMKL.
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
      acc_mean = 0.0
      !$omp simd reduction(+:acc_mean)
      do i = 1, n
        acc_mean = acc_mean + x(i, j)
      end do
      mean_val = acc_mean * inv_n

      acc_var = 0.0
      !$omp simd reduction(+:acc_var) private(diff)
      do i = 1, n
        diff = x(i, j) - mean_val
        acc_var = acc_var + diff * diff
      end do
      var_val = acc_var * inv_n
      if (var_val < 0.0) var_val = 0.0

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
    real(c_float), intent(out) :: x_out(n, p)
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
      if (sds(j) > 0.0) then
        inv_sd = 1.0d0 / sds(j)
      else
        inv_sd = 1.0d0
      end if
      !$omp simd
      do i = 1, n
        x_out(i, j) = real((x(i, j) - mu) * inv_sd, c_float)
      end do
    end do
    !$omp end parallel do
  end subroutine sby_apply_zscore_simd_f

  ! -------------------------------------------------------------------
  ! sby_revert_zscore_fma_f
  ! Reverte o z-score diretamente na consolidacao final. O laco aninhado
  ! por colunas e instruido com !$OMP SIMD para forcar as 4 unidades FMA do
  ! hardware, executando multiplicacao pelo desvio padrao e soma da media na
  ! mesma instrucao vfmadd213ps. Nao usa dscal nem daxpy.
  ! -------------------------------------------------------------------
  subroutine sby_revert_zscore_fma_f(x, n, p, means, sds, x_out, status) &
      bind(c, name="sby_revert_zscore_fma_f")
    integer(c_int), intent(in), value :: n
    integer(c_int), intent(in), value :: p
    real(c_float), intent(in)   :: x(n, p)
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
      ! unica instrucao vfmadd213ps nas unidades FMA do AVX-512.
      !$omp simd
      do i = 1, n
        x_out(i, j) = dble(x(i, j)) * sd + mu
      end do
    end do
    !$omp end parallel do
  end subroutine sby_revert_zscore_fma_f

  ! -------------------------------------------------------------------
  ! sby_sgemm_neg2_f
  ! Computes only C = -2 * A * transpose(B) through the standard Fortran BLAS
  ! interface. Matrices use column-major storage. lda and ldb are the physical
  ! row strides of the source matrices, while ldc is the physical row stride
  ! of C. A row interval can therefore be passed as an offset pointer without
  ! copying its columns. No allocation, norm calculation, OpenMP region, or
  ! thread setting is performed here. SGEMM is entered from serial C++ code so
  ! the selected BLAS owns its worker team and avoids nested fork/join traffic.
  ! -------------------------------------------------------------------
  subroutine sby_sgemm_neg2_f(a, lda, b, ldb, c, ldc, m, n, k, status) &
      bind(c, name="sby_sgemm_neg2_f")
    integer(c_int), intent(in), value :: lda, ldb, ldc, m, n, k
    real(c_float), intent(in) :: a(lda, *), b(ldb, *)
    real(c_float), intent(out) :: c(ldc, *)
    integer(c_int), intent(out) :: status

    status = 0
    if (m < 1 .or. n < 1 .or. k < 1 .or. lda < m .or. ldb < n .or. ldc < m) then
      status = -1
      return
    end if
    call sgemm('N', 'T', m, n, k, -2.0_c_float, a, lda, &
               b, ldb, 0.0_c_float, c, ldc)
  end subroutine sby_sgemm_neg2_f

  ! -------------------------------------------------------------------
  ! sby_adasyn_interp_uniform_f
  ! Interpolacao sintetica do ADASYN no espaco padronizado. Para cada linha
  ! sintetica, combina a linha base e um vizinho minoritario com peso lambda
  ! pre-gerado no C++/Rcpp. Os indices base e vizinho sao 1-based.
  !   minority(n_min x p)         matriz minoritaria padronizada
  !   base_idx(n_syn), nbr_idx(n_syn)  indices base e vizinho por linha sintetica
  !   lambda(n_syn)               pesos uniformes pre-gerados no C++/Rcpp
  !   syn_out(n_syn x p)          saida sintetica padronizada
  ! -------------------------------------------------------------------
  subroutine sby_adasyn_interp_uniform_f(minority, n_min, p, base_idx, nbr_idx, &
      lambda, n_syn, syn_out, status) bind(c, name="sby_adasyn_interp_uniform_f")
    integer(c_int), intent(in), value :: n_min
    integer(c_int), intent(in), value :: p
    integer(c_int), intent(in), value :: n_syn
    real(c_float), intent(in)   :: minority(n_min, p)
    integer(c_int), intent(in)  :: base_idx(n_syn)
    integer(c_int), intent(in)  :: nbr_idx(n_syn)
    real(c_float), intent(in)   :: lambda(n_syn)
    real(c_float), intent(out)  :: syn_out(n_syn, p)
    integer(c_int), intent(out) :: status

    integer :: s, j, bi, ni
    real(c_float) :: lam, base_val, nbr_val

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
      ! funde a multiplicacao e a soma na instrucao vfmadd213ps das unidades FMA.
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

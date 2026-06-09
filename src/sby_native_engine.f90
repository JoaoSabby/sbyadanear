module sby_native_engine_mkl_mod
  use, intrinsic :: iso_c_binding
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

  subroutine sby_compute_zscore_population_f(x, n, p, means, sds, status) bind(c, name="sby_compute_zscore_population_f")
    integer(c_int), intent(in), value :: n
    integer(c_int), intent(in), value :: p
    real(c_double), intent(in)  :: x(n, p)
    real(c_double), intent(out) :: means(p)
    real(c_double), intent(out) :: sds(p)
    integer(c_int), intent(out) :: status

    integer :: j
    real(c_double) :: inv_n, mean_val, ex2, sumsq

    status = 0
    if (n < 1 .or. p < 1) then
      status = -1
      return
    end if

    inv_n = 1.0d0 / dble(n)

    ! === 1. Soma de medias por coluna contigua do layout R n x p ===
    !$omp parallel do default(none) shared(x, means, p, n, inv_n) &
    !$omp& private(j, mean_val) schedule(static)
    do j = 1, p
      mean_val = sum(x(:, j))
      means(j) = mean_val * inv_n
    end do
    !$omp end parallel do

    ! === 2. Soma dos quadrados via CBLAS com stride 1 em colunas contiguas ===
    !$omp parallel do default(none) shared(x, sds, means, p, n, inv_n) &
    !$omp& private(j, ex2, mean_val, sumsq) schedule(static)
    do j = 1, p
      ex2 = cblas_ddot(n, x(1, j), 1, x(1, j), 1)
      mean_val = means(j)
      sumsq = ex2 * inv_n - mean_val * mean_val
      if (sumsq < 0.0d0 .and. sumsq > -1.0d-12) sumsq = 0.0d0
      if (sumsq < 0.0d0) sumsq = 0.0d0
      sds(j) = sqrt(sumsq)
    end do
    !$omp end parallel do
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

    status = 0
    if (n < 1 .or. p < 1) then
      status = -1
      return
    end if

    !$omp parallel do default(none) shared(x, x_out, means, sds, p, n) &
    !$omp& private(j, i, alpha, mu) schedule(static)
    do j = 1, p
      call cblas_dcopy(n, x(1, j), 1, x_out(1, j), 1)
      mu = means(j)
      do i = 1, n
        x_out(i, j) = x_out(i, j) - mu
      end do
      if (sds(j) > 0.0d0) then
        alpha = 1.0d0 / sds(j)
      else
        alpha = 1.0d0
      end if
      call cblas_dscal(n, alpha, x_out(1, j), 1)
    end do
    !$omp end parallel do
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

    status = 0
    if (n < 1 .or. p < 1) then
      status = -1
      return
    end if

    !$omp parallel do default(none) shared(x, x_out, means, sds, p, n) &
    !$omp& private(j, i, alpha, mu) schedule(static)
    do j = 1, p
      call cblas_dcopy(n, x(1, j), 1, x_out(1, j), 1)
      alpha = sds(j)
      call cblas_dscal(n, alpha, x_out(1, j), 1)
      mu = means(j)
      do i = 1, n
        x_out(i, j) = x_out(i, j) + mu
      end do
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

    status = 0
    if (p < 1 .or. n1 < 0 .or. n2 < 0) then
      status = -1
      return
    end if

    !$omp parallel do default(none) shared(a, b, c_out, n1, n2, p) private(j) schedule(static)
    do j = 1, p
      if (n1 > 0) then
        call cblas_dcopy(n1, a(1, j), 1, c_out(1, j), 1)
      end if
      if (n2 > 0) then
        call cblas_dcopy(n2, b(1, j), 1, c_out(n1 + 1, j), 1)
      end if
    end do
    !$omp end parallel do
  end subroutine sby_rbind_matrix_f

end module sby_native_engine_mkl_mod

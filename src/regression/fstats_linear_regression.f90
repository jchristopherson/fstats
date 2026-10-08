module fstats_linear_regression
    use iso_fortran_env
    use fstats_regression
    use blas
    use fstats_errors
    use fstats_distributions
    implicit none
    private
    public :: linear_least_squares

contains
! ------------------------------------------------------------------------------
subroutine linear_least_squares(order, intercept, x, y, coeffs, &
    ymod, resid, stats, alpha)
    !! Computes a linear least-squares regression to fit a set of data.
    !!
    !! See Also
    !!
    !! - <a href="https://en.wikipedia.org/wiki/Linear_regression" target="_blank">Wikipedia - Linear Regression</a>
    !! - <a href="https://www.spcforexcel.com/knowledge/root-cause-analysis/understanding-regression-statistics-part-1" 
    !! target="_blank">SPC Excel Understanding Regression Statistics</a>
    integer(int32), intent(in) :: order
        !! The order of the equation to fit.  This value must be at 
        !! least one (linear equation), but can be higher as desired, 
        !! as long as there is sufficient data.
    logical, intent(in) :: intercept
        !! Set to true if the intercept is being computed as part of 
        !! the regression; else, false.
    real(real64), intent(in) :: x(:)
        !! An N-element array containing the independent variable
        !! measurement points.
    real(real64), intent(in) :: y(:)
        !! An N-element array containing the dependent variable
        !! measurement points.
    real(real64), intent(out) :: coeffs(:)
        !! An ORDER+1 element array where the coefficients will be written.
    real(real64), intent(out) :: ymod(:)
        !! An N-element array where the modeled data will be written.
    real(real64), intent(out) :: resid(:)
        !! An N-element array where the residual error data will be 
        !! written (modeled - actual).
    type(regression_statistics), intent(out), optional :: stats(:)
        !! An M-element array of regression_statistics items where 
        !! M = ORDER + 1 when intercept is set to true; however, if 
        !! intercept is set to false, M = ORDER.
    real(real64), intent(in), optional :: alpha
        !! The significance level at which to evaluate the confidence 
        !! intervals.  The default value is 0.05 such that a 95% 
        !! confidence interval is calculated.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: half = 0.5d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, npts, ncols, ncoeffs
    real(real64) :: alph, var, df, ssr, talpha
    real(real64), allocatable :: a(:,:), c(:,:), cxt(:,:)
    type(t_distribution) :: dist
    
    ! Initialization
    npts = size(x)
    ncoeffs = order + 1
    ncols = order
    if (intercept) ncols = ncols + 1
    alph = 0.05d0
    if (present(alpha)) alph = alpha

    ! Input Check
    if (order < 1) error stop FS_INVALID_INPUT_ERROR
    if (size(y) /= npts) error stop FS_ARRAY_SIZE_ERROR
    if (size(coeffs) /= ncoeffs) error stop FS_ARRAY_SIZE_ERROR
    if (size(ymod) /= npts) error stop FS_ARRAY_SIZE_ERROR
    if (size(resid) /= npts) error stop FS_ARRAY_SIZE_ERROR
    if (present(stats)) then
        if (size(stats) /= ncols) error stop FS_ARRAY_SIZE_ERROR
    end if

    ! Memory Allocation
    allocate(a(npts, ncols), c(ncols, ncols), cxt(ncols, npts))

    ! Compute the coefficient matrix
    a = design_matrix(order, intercept, x)

    ! Compute the covariance matrix
    c = covariance_matrix(a)

    ! Compute the coefficients (NCOLS-by-1)
    call DGEMM("N", "T", ncols, npts, ncols, one, c, ncols, a, npts, zero, &
        cxt, ncols)     ! C * X**T

    i = 2
    coeffs(1) = zero
    if (intercept) i = 1
    call DGEMM("N", "N", ncols, 1, npts, one, cxt, ncols, y, npts, zero, &
        coeffs(i:), ncols)  ! (C * X**T) * Y

    ! Evaluate the model and compute the residuals
    call DGEMM("N", "N", npts, 1, ncols, one, a, npts, coeffs(i:), &
        ncols, zero, ymod, npts)
    resid = ymod - y

    ! If the user doesn't want the statistics calculations we can stop now
    if (.not.present(stats)) return
    
    ! Start the process of computing statistics
    stats = calculate_regression_statistics(resid, coeffs(i:), c, alph)
end subroutine

! ------------------------------------------------------------------------------
end module
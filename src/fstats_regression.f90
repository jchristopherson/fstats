module fstats_regression
    use iso_fortran_env
    use linalg
    use fstats_errors
    use blas
    use ferror
    use fstats_descriptive_statistics
    use fstats_distributions
    use fstats_special_functions
    use fstats_hypothesis
    implicit none
    private
    public :: iteration_controls
    public :: convergence_info
    public :: lm_solver_options
    public :: regression_function
    public :: iteration_update
    public :: regression_statistics
    public :: r_squared
    public :: adjusted_r_squared
    public :: correlation
    public :: design_matrix
    public :: covariance_matrix
    public :: linear_least_squares
    public :: calculate_regression_statistics
    public :: jacobian
    public :: nonlinear_least_squares
    public :: FS_LEVENBERG_MARQUARDT_UPDATE
    public :: FS_QUADRATIC_UPDATE
    public :: FS_NIELSEN_UPDATE

! ******************************************************************************
! CONSTANTS
! ------------------------------------------------------------------------------
    integer(int32), parameter :: FS_LEVENBERG_MARQUARDT_UPDATE = 1
    integer(int32), parameter :: FS_QUADRATIC_UPDATE = 2
    integer(int32), parameter :: FS_NIELSEN_UPDATE = 3

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    type regression_statistics
       !! A container for regression-related statistical information. 
        real(real64) :: standard_error
            !! The standard error for the model coefficient.
            !!
            !! $$ E_{s}(\beta_{i}) = \sqrt{\sigma^{2} C_{ii}} $$
        real(real64) :: t_statistic
            !! The T-statistic for the model coefficient.
            !!
            !! $$ t_o = \frac{ \beta_{i} }{E_{s}(\beta_{i})} $$
        real(real64) :: probability
            !! The probability that the coefficient is not statistically 
            !! important.  A statistically important coefficient will have a 
            !! low probability (p-value), typically 0.05 or lower; however, a 
            !! p-value of up to ~0.2 may be acceptable dependent upon the 
            !! problem.  Typically any p-value larger than ~0.2 indicates the 
            !! parameter is not statistically important for the model.
            !!
            !! $$ p = t_{|t_o|, df_{residual}} $$
        real(real64) :: confidence_interval
            !! The confidence interval for the parameter at the level 
            !! determined by the regression process.
            !!
            !! $$ c = t_{\alpha, df} E_{s}(\beta_{i}) $$
    end type

    type iteration_controls
        !! Provides a collection of iteration control parameters.
        integer(int32) :: max_iteration_count
            !! Defines the maximum number of iterations allowed.
        integer(int32) :: max_function_evaluations
            !! Defines the maximum number of function evaluations allowed.
        real(real64) :: gradient_tolerance
            !! Defines a tolerance on the gradient of the fitted function.
        real(real64) :: change_in_solution_tolerance
            !! Defines a tolerance on the change in parameter values.
        real(real64) :: residual_tolerance
            !! Defines a tolerance on the metric associated with the residual 
            !! error.
        real(real64) :: iteration_improvement_tolerance
            !! Defines a tolerance to ensure adequate improvement on each 
            !! iteration.
        integer(int32) :: max_iteration_between_updates
            !! Defines how many iterations can pass before a re-evaluation of 
            !! the Jacobian matrix is forced.
    contains
        procedure, public :: set_to_default => lm_set_default_tolerances
    end type

    type convergence_info
        !! Provides information regarding convergence status.
        logical :: converge_on_gradient
            !! True if convergence on the gradient was achieved; else, false.
        real(real64) :: gradient_value
            !! The value of the gradient test parameter.
        logical :: converge_on_solution_change
            !! True if convergence on the change in solution was achieved; else,
            !! false.
        real(real64) :: solution_change_value
            !! The value of the change in solution parameter.
        logical :: converge_on_residual_parameter
            !! True if convergence on the residual error parameter was achieved; 
            !! else, false.
        real(real64) :: residual_value
            !! The value of the residual error parameter.
        logical :: reach_iteration_limit
            !! True if the solution did not converge in the allowed number of 
            !! iterations.
        integer(int32) :: iteration_count
            !! The iteration count.
        logical :: reach_function_evaluation_limit
            !! True if the solution did not converge in the allowed number of
            !! function evaluations.
        integer(int32) :: function_evaluation_count
            !! The function evaluation count.
        logical :: user_requested_stop
            !! True if the user requested the stop; else, false.
    end type

    type lm_solver_options
        !! Options to control the Levenberg-Marquardt solver.
        integer(int32) :: method
            !! The solver method to utilize.
            !! - FS_LEVENBERG_MARQUARDT_UPDATE:
            !! - FS_QUADRATIC_UPDATE:
            !! - FS_NIELSEN_UDPATE:
        real(real64) :: finite_difference_step_size
            !! The step size used for the finite difference calculations of the
            !! Jacobian matrix.
        real(real64) :: damping_increase_factor
            !! The factor to use when increasing the damping parameter.
        real(real64) :: damping_decrease_factor
            !! The factor to use when decreasing the damping parameter.
    contains
        procedure, public :: set_to_default => lm_set_default_settings
    end type

    interface
        subroutine regression_function(xdata, params, f, stop)
            !! Defines the interface of a subroutine computing the function
            !! values at each of the N data points as part of a regression
            !! analysis.
            use iso_fortran_env, only : real64
            real(real64), intent(in), dimension(:) :: xdata
                !! An N-element array containing the N independent data points.
            real(real64), intent(in), dimension(:) :: params
                !! An M-element array containing the M model parameters.
            real(real64), intent(out), dimension(:) :: f
                !! An N-element array where the results of the N function 
                !! evaluations will be written.
            logical, intent(out) :: stop
                !! A mechanism to force a stop to the iteration process.  If
                !! set to true, the iteration process will terminate.  If set
                !! to false, the iteration process will continue along as 
                !! normal.
        end subroutine

        subroutine iteration_update(iter, funvals, resid, params, step)
            !! Defines a routine for providing updates about an iteration
            !! process.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: iter
                !! The current iteration number.
            real(real64), intent(in), dimension(:) :: funvals
                !! The function values.
            real(real64), intent(in), dimension(:) :: resid
                !! The residuals.
            real(real64), intent(in), dimension(:) :: params
                !! The model parameters.
            real(real64), intent(in), dimension(:) :: step
                !! Step sizes for each parameter.
        end subroutine
    end interface

contains

! ------------------------------------------------------------------------------
function r_squared(x, xm, err) result(rst)
    !! Computes the R-squared value for a data set.
    !!
    !! The R-squared value is computed by determining the sum of the squares
    !! of the residuals: 
    !! $$ SS_{res} = \Sigma \left( y_i - f_i \right)^2 $$
    !! The total sum of the squares: 
    !! $$ SS_{tot} = \Sigma \left( y_i - \bar{y} \right)^2 $$. 
    !! The R-squared value is then: 
    !! $$ R^2 = 1 - \frac{SS_{res}}{SS_{tot}} $$.
    !!
    !! See Also:
    !!
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Coefficient_of_determination)
    real(real64), intent(in) :: x(:)
        !! An N-element array containing the dependent variables from 
        !! the data set.
    real(real64), intent(in) :: xm(:)
        !! An N-element array containing the corresponding modeled 
        !! values.
    class(errors), intent(inout), optional, target :: err
        !! A mechanism for communicating errors and warnings
        !! to the caller.  Possible warning and error codes are as 
        !! follows.
        !! - FS_NO_ERROR: No errors encountered.
        !! - FS_ARRAY_SIZE_ERROR: Occurs if x and xm are not the 
        !!   same size.
    real(real64) :: rst
        !! The result.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: esum, vt
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    n = size(x)
    if (size(xm) /= n) then
        call report_array_size_error(errmgr, "r_squared_real64", "XM", n, &
            size(xm))
        return
    end if

    ! Process
    esum = zero
    do i = 1, n
        esum = esum + (x(i) - xm(i))**2
    end do
    vt = variance(x) * (n - one)
    rst = one - esum / vt
end function

! ------------------------------------------------------------------------------
function adjusted_r_squared(p, x, xm, err) result(rst)
    !! Computes the adjusted R-squared value for a data set.
    !!
    !! The adjusted R-squared provides a mechanism for tempering the effects
    !! of extra explanatory variables on the traditional R-squared 
    !! calculation.  It is computed by noting the sample size \( n \) and 
    !! the number of variables \( p \).
    !! $$ \bar{R}^2 = 1 - \left( 1 - R^2 \right) \frac{n - 1}{n - p} $$.
    !!
    !! See Also:
    !!
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Coefficient_of_determination#Adjusted_R2)
    integer(int32), intent(in) :: p
        !! The number of variables.
    real(real64), intent(in) :: x(:)
        !! An N-element array containing the dependent variables from 
        !! the data set.
    real(real64), intent(in) :: xm(:)
        !! An N-element array containing the corresponding modeled 
        !! values.
    class(errors), intent(inout), optional, target :: err
        !! A mechanism for communicating errors and warnings
        !! to the caller.  Possible warning and error codes are as 
        !! follows.
        !! - FS_NO_ERROR: No errors encountered.
        !! - FS_ARRAY_SIZE_ERROR: Occurs if x and xm are not the 
        !!   same size.
    real(real64) :: rst
        !! The result.

    ! Local Variables
    integer(int32) :: n
    real(real64) :: r2
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Parameters
    real(real64), parameter :: one = 1.0d0
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)

    ! Process
    r2 = r_squared(x, xm, errmgr)
    if (errmgr%has_error_occurred()) return
    rst = one - (one - r2) * (n - one) / (n - p - one)
end function

! ------------------------------------------------------------------------------
pure function correlation(x, y) result(rst)
    !! Computes the sample correlation coefficient (an estimate to the 
    !! population Pearson correlation) as follows.
    !!
    !! $$ r_{xy} = \frac{cov(x, y)}{s_{x} s_{y}} $$.
    !!
    !! Where, \( s_{x} \) & \( s_{y} \) are the sample standard deviations of
    !! x and y respectively.
    real(real64), intent(in), dimension(:) :: x
        !! The first N-element data set.
    real(real64), intent(in), dimension(size(x)) :: y
        !! The second N-element data set.
    real(real64) :: rst
        !! The correlation coefficient.

    ! Process
    rst = covariance(x, y) / (standard_deviation(x) * standard_deviation(y))
end function

! ------------------------------------------------------------------------------
subroutine design_matrix(order, intercept, x, c, err)
    !! Computes the design matrix \( X \) for the linear 
    !! least-squares regression problem of \( X \beta = y \), where 
    !! \( X \) is the matrix computed here, \( \beta \) is 
    !! the vector of coefficients to be determined, and \( y \) is the 
    !! vector of measured dependent variables.
    !!
    !! See Also
    !!
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Linear_regression)
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Vandermonde_matrix)
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Design_matrix)
    integer(int32), intent(in) :: order
        !! The order of the equation to fit.  This value must be
        !! at least one (linear equation), but can be higher as desired.
    logical, intent(in) :: intercept
        !! Set to true if the intercept is being computed
        !! as part of the regression; else, false.
    real(real64), intent(in) :: x(:)
        !! An N-element array containing the independent variable
        !! measurement points.
    real(real64), intent(out) :: c(:,:)
        !! An N-by-K matrix where the results will be written.  K
        !! must equal order + 1 in the event intercept is true; 
        !! however, if intercept is false, K must equal order.
    class(errors), intent(inout), optional, target :: err
        !! A mechanism for communicating errors and warnings to the 
        !! caller.  Possible warning and error codes are as follows.
        !! - FS_NO_ERROR: No errors encountered.
        !! - FS_ARRAY_SIZE_ERROR: Occurs if c is not properly sized.
        !! - FS_INVALID_INPUT_ERROR: Occurs if order is less than 1.

    ! Parameters
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, start, npts, ncols
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    npts = size(x)
    ncols = order
    if (intercept) ncols = ncols + 1

    ! Input Check
    if (order < 1) then
        call errmgr%report_error("design_matrix", &
            "The model order must be at least one.", FS_INVALID_INPUT_ERROR)
        return
    end if
    if (size(c, 1) /= npts .or. size(c, 2) /= ncols) then
        call report_matrix_size_error(errmgr, "design_matrix", &
            "c", npts, ncols, size(c, 1), size(c, 2))
        return
    end if

    ! Process
    if (intercept) then
        c(:,1) = one
        c(:,2) = x
        start = 3
    else
        c(:,1) = x
        start = 2
    end if
    if (start >= ncols) return
    do i = start, ncols
        c(:,i) = c(:,i-1) * x
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine covariance_matrix(x, c, err)
    !! Computes the covariance matrix \( C \) where 
    !! \( C = \left( X^{T} X \right)^{-1} \) and \( X \) is computed
    !! by design_matrix.
    !!
    !! See Also
    !!
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Covariance_matrix)
    !! - [Wikipedia - Regression](https://en.wikipedia.org/wiki/Linear_regression)
    real(real64), intent(in) :: x(:,:)
        !! An M-by-N matrix containing the formatted independent data
        !!  matrix \( X \) as computed by design_matrix.
    real(real64), intent(out) :: c(:,:)
        !! The N-by-N covariance matrix.
    class(errors), intent(inout), optional, target :: err
        !! A mechanism for communicating errors and warnings to the 
        !! caller.  Possible warning and error codes are as follows.
        !! - FS_NO_ERROR: No errors encountered.
        !! - FS_ARRAY_SIZE_ERROR: Occurs if any of the matrices are not 
        !!      sized correctly.
        !! - FS_MEMORY_ERROR: Occurs if there is a memory allocation 
        !!      error.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    integer(int32) :: npts, ncoeffs, flag
    real(real64), allocatable :: xtx(:,:)
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    npts = size(x, 1)
    ncoeffs = size(x, 2)

    ! Input Checking
    if (size(c, 1) /= ncoeffs .or. size(c, 2) /= ncoeffs) then
        call report_matrix_size_error(errmgr, "covariance_matrix", &
            "c", ncoeffs, ncoeffs, size(c, 1), size(c, 2))
        return
    end if

    ! Local Memory Allocation
    allocate(xtx(ncoeffs, ncoeffs), stat = flag)
    if (flag /= 0) then
        call report_memory_error(errmgr, "covariance_matrix", flag)
        return
    end if

    ! Compute X**T * X
    call DGEMM("T", "N", ncoeffs, ncoeffs, npts, one, x, npts, x, npts, &
        zero, xtx, ncoeffs)
    
    ! Compute the inverse of X**T * X to obtain the covariance matrix
    call mtx_pinverse(xtx, c, err = errmgr)
    if (errmgr%has_error_occurred()) return
end subroutine

! ------------------------------------------------------------------------------
subroutine linear_least_squares(order, intercept, x, y, coeffs, &
    ymod, resid, stats, alpha, err)
    !! Computes a linear least-squares regression to fit a set of data.
    !!
    !! See Also
    !!
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Linear_regression)
    !! - [SPC Excel Understanding Regression Statistics](https://www.spcforexcel.com/knowledge/root-cause-analysis/understanding-regression-statistics-part-1)
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
    class(errors), intent(inout), optional, target :: err
        !! A mechanism for communicating errors and warnings to the 
        !! caller.  Possible warning and error codes are as follows.
        !! - FS_NO_ERROR: No errors encountered.
        !! - FS_ARRAY_SIZE_ERROR: Occurs if any of the arrays are not 
        !!      approriately sized.
        !! - FS_INVALID_INPUT_ERROR: Occurs if order is less than 1.
        !! - FS_MEMORY_ERROR: Occurs if there is a memory allocation 
        !!      error.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: half = 0.5d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, npts, ncols, ncoeffs, flag
    real(real64) :: alph, var, df, ssr, talpha
    real(real64), allocatable :: a(:,:), c(:,:), cxt(:,:)
    type(t_distribution) :: dist
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    npts = size(x)
    ncoeffs = order + 1
    ncols = order
    if (intercept) ncols = ncols + 1
    alph = 0.05d0
    if (present(alpha)) alph = alpha

    ! Input Check
    if (order < 1) then
        call errmgr%report_error("linear_least_squares", &
            "The model order must be at least one.", FS_INVALID_INPUT_ERROR)
        return
    end if
    if (size(y) /= npts) then
        call report_array_size_error(errmgr, "linear_least_squares", &
            "y", npts, size(y))
        return
    end if
    if (size(coeffs) /= ncoeffs) then
        call report_array_size_error(errmgr, "linear_least_squares", &
            "coeffs", ncoeffs, size(coeffs))
        return
    end if
    if (size(ymod) /= npts) then
        call report_array_size_error(errmgr, "linear_least_squares", &
            "ymod", npts, size(ymod))
        return
    end if
    if (size(resid) /= npts) then
        call report_array_size_error(errmgr, "linear_least_squares", &
            "resid", npts, size(resid))
        return
    end if
    if (present(stats)) then
        if (size(stats) /= ncols) then
            call report_array_size_error(errmgr, &
                "linear_least_squares", "stats", ncols, size(stats))
            return
        end if
    end if

    ! Memory Allocation
    allocate(a(npts, ncols), stat = flag)
    if (flag == 0) allocate(c(ncols, ncols), stat = flag)
    if (flag == 0) allocate(cxt(ncols, npts), stat = flag)
    if (flag /= 0) then
        call report_memory_error(errmgr, "linear_least_squares", flag)
        return
    end if

    ! Compute the coefficient matrix
    call design_matrix(order, intercept, x, a, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute the covariance matrix
    call covariance_matrix(a, c, errmgr)
    if (errmgr%has_error_occurred()) return

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
    stats = calculate_regression_statistics(resid, coeffs(i:), c, alph, &
        errmgr)
end subroutine

! ------------------------------------------------------------------------------
function calculate_regression_statistics(resid, params, c, alpha, err) &
    result(rst)
    !! Computes statistics for the quality of fit for a regression 
    !! model.
    real(real64), intent(in) :: resid(:)
        !! An M-element array containing the model residual errors.
    real(real64), intent(in) :: params(:)
        !! An N-element array containing the model parameters.
    real(real64), intent(in) :: c(:,:)
        !! The N-by-N covariance matrix.
    real(real64), intent(in), optional :: alpha
        !! The significance level at which to evaluate the confidence 
        !! intervals.  The default value is 0.05 such that a 95% 
        !! confidence interval is calculated.
    class(errors), intent(inout), optional, target :: err
        !! A mechanism for communicating errors and warnings to the 
        !! caller.  Possible warning and error codes are as follows.
        !! - FS_NO_ERROR: No errors encountered.
        !! - FS_ARRAY_SIZE_ERROR: Occurs if c is not sized correctly.
        !! - FS_INVALID_INPUT_ERROR: Occurs if order is less than 1.
        !! - FS_MEMORY_ERROR: Occurs if there is a memory allocation 
        !!      error.
    type(regression_statistics), allocatable :: rst(:)
        !! A regression_statistics object containing the analysis results.

    ! Parameters
    real(real64), parameter :: p05 = 0.05d0
    real(real64), parameter :: half = 0.5d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, m, n, dof, flag
    real(real64) :: a, ssr, var, talpha
    type(t_distribution) :: dist
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Initialization
    m = size(resid)
    n = size(params)
    dof = m - n
    if (present(alpha)) then
        a = alpha
    else
        a = p05
    end if
    allocate(rst(n), stat = flag)
    if (flag /= 0) then
    end if

    ! Input Checking
    if (size(c, 1) /= n .or. size(c, 2) /= n) then
    end if

    ! Process
    ssr = norm2(resid)**2   ! sum of the squares of the residual
    var = ssr / dof
    dist%dof = real(dof, real64)
    talpha = confidence_interval(dist, a, one, 1)
    do i = 1, n
        rst(i)%standard_error = sqrt(var * c(i,i))
        rst(i)%t_statistic = params(i) / rst(i)%standard_error
        rst(i)%probability = regularized_beta( &
            half * dof, &
            half, &
            real(dof, real64) / (dof + (rst(i)%t_statistic)**2) &
        )
        rst(i)%confidence_interval = talpha * rst(i)%standard_error
    end do
end function

! ------------------------------------------------------------------------------
subroutine jacobian(fun, xdata, params, &
    jac, stop, f0, f1, step, err)
    !! Computes the Jacobian matrix for a nonlinear regression problem.
    procedure(regression_function), intent(in), pointer :: fun
        !! A pointer to the regression_function to evaluate.
    real(real64), intent(in) :: xdata(:)
        !! The M-element array containing x-coordinate data.
    real(real64), intent(in) :: params(:)
        !! The N-element array containing the model parameters.
    real(real64), intent(out) :: jac(:,:)
        !! The M-by-N matrix where the Jacobian will be written.
    logical, intent(out) :: stop
        !! A value that the user can set in fun forcing the
        !! evaluation process to stop prior to completion.
    real(real64), intent(in), optional, target :: f0(:)
        !! An optional M-element array containing the model values
        !!  using the current parameters as defined in m.  This input 
        !! can be used to prevent the routine from performing a 
        !! function evaluation at the model parameter state defined in 
        !! params.
    real(real64), intent(out), optional, target :: f1(:)
        !! An optional M-element workspace array used for function
        !! evaluations.
    real(real64), intent(in), optional :: step
        !! The differentiation step size.  The default is the square 
        !! root of machine precision.
    class(errors), intent(inout), optional, target :: err
        !! A mechanism for communicating errors and warnings to the 
        !! caller.  Possible warning and error codes are as follows.
        !! - FS_NO_ERROR: No errors encountered.
        !! - FS_ARRAY_SIZE_ERROR: Occurs if any of the arrays are not 
        !!      properly sized.
        !! - FS_MEMORY_ERROR: Occurs if there is a memory allocation 
        !!      error.

    ! Local Variables
    real(real64) :: h
    integer(int32) :: m, n, flag, expected, actual
    real(real64), pointer :: f1p(:), f0p(:)
    real(real64), allocatable, target :: f1a(:), f0a(:), work(:)
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (present(step)) then
        h = step
    else
        h = sqrt(epsilon(h))
    end if
    m = size(xdata)
    n = size(params)

    ! Input Size Checking
    if (size(jac, 1) /= m .or. size(jac, 2) /= n) then
        call report_matrix_size_error(errmgr, "jacobian", &
            "JAC", m, n, size(jac, 1), size(jac, 2))
        return
    end if
    if (present(f0)) then
        ! Check Size
        if (size(f0) /= m) then
            call report_array_size_error(errmgr, "jacobian", &
                "F0", m, size(f0))
            return
        end if
        f0p(1:m) => f0
    else
        ! Allocate space, and fill the array with the current function
        ! results
        allocate(f0a(m), stat = flag)
        if (flag /= 0) go to 20
        f0p(1:m) => f0a
        call fun(xdata, params, f0p, stop)
        if (stop) return
    end if
    if (present(f1)) then
        ! Check Size
        if (size(f1) /= m) then
            call report_array_size_error(errmgr, "jacobian", &
                "F1", m, size(f1))
            return
        end if
        f1p(1:m) => f1
    else
        ! Allocate space
        allocate(f1a(m), stat = flag)
        if (flag /= 0) go to 20
        f1p(1:m) => f1a
    end if

    ! Allocate a workspace array the same size as params
    allocate(work(n), stat = flag)
    if (flag /= 0) go to 20

    ! Compute the Jacobian
    call jacobian_finite_diff(fun, xdata, params, f0p, jac, f1p, &
        stop, h, work)

    ! End
    return

    ! Memroy Allocation Error Handling
20  continue
    call report_memory_error(errmgr, "jacobian", flag)
    return
end subroutine

! ------------------------------------------------------------------------------
subroutine nonlinear_least_squares(fun, x, y, params, ymod, &
    resid, weights, maxp, minp, stats, alpha, controls, settings, info, &
    status, cov, err)
    !! Performs a nonlinear regression to fit a model using a version
    !! of the Levenberg-Marquardt algorithm.
    procedure(regression_function), intent(in), pointer :: fun
        !! A pointer to the regression_function to evaluate.
    real(real64), intent(in) :: x(:)
        !! The M-element array containing independent data.
    real(real64), intent(in) :: y(:)
        !! The M-element array containing dependent data.
    real(real64), intent(inout) :: params(:)
        !! On input, the N-element array containing the initial estimate
        !! of the model parameters.  On output, the computed model 
        !! parameters.
    real(real64), intent(out) :: ymod(:)
        !! An M-element array where the modeled dependent data will
        !! be written.
    real(real64), intent(out) :: resid(:)
        !! An M-element array where the model residuals will be
        !! written.
    real(real64), intent(in), optional, target :: weights(:)
        !! An optional M-element array allowing the weighting of
        !! individual points.
    real(real64), intent(in), optional, target :: maxp(:)
        !! An optional N-element array that can be used as upper limits 
        !! on the parameter values.  If no upper limit is requested for
        !! a particular parameter, utilize a very large value.  The 
        !! internal default is to utilize huge() as a value.
    real(real64), intent(in), optional, target :: minp(:)
        !! An optional N-element array that can be used as lower limits 
        !! on the parameter values.  If no lower limit is requested for
        !! a particalar parameter, utilize a very large magnitude, but 
        !! negative, value.  The internal default is to utilize -huge() 
        !! as a value.
    type(regression_statistics), intent(out), optional :: stats(:)
        !! An optional N-element array that, if supplied, will be used 
        !! to return statistics about the fit for each parameter.
    real(real64), intent(in), optional :: alpha
        !! The significance level at which to evaluate the confidence 
        !! intervals.  The default value is 0.05 such that a 95% 
        !! confidence interval is calculated.
    type(iteration_controls), intent(in), optional :: controls
        !! An optional input providing custom iteration controls.
    type(lm_solver_options), intent(in), optional :: settings
        !! An optional input providing custom settings for the solver.
    type(convergence_info), intent(out), optional, target :: info
        !! An optional output that can be used to gain information about
        !! the iterative solution and the nature of the convergence.
    procedure(iteration_update), intent(in), pointer, optional :: status
        !! An optional pointer to a routine that can be used to extract
        !! iteration information.
    real(real64), intent(out), optional, dimension(:,:) :: cov
        !! An optional N-by-N matrix that, if supplied, will be used to return
        !! the covariance matrix.
    class(errors), intent(inout), optional, target :: err
        !! A mechanism for communicating errors and warnings to the 
        !! caller.  Possible warning and error codes are as follows.
        !!
        !! - FS_NO_ERROR: No errors encountered.
        !! - FS_ARRAY_SIZE_ERROR: Occurs if any of the arrays are not 
        !!      properly sized.
        !!
        !! - FS_MEMORY_ERROR: Occurs if there is a memory allocation 
        !!      error.
        !!
        !! - FS_UNDERDEFINED_PROBLEM_ERROR: Occurs if the problem posed 
        !!      is underdetetermined (M < N).
        !!
        !! - FS_TOLERANCE_TOO_SMALL_ERROR: Occurs if any supplied 
        !!      tolerances are too small to be practical.
        !!
        !! - FS_TOO_FEW_ITERATION_ERROR: Occurs if too few iterations 
        !!      are allowed.

    ! Parameters
    real(real64), parameter :: too_small = 1.0d-14
    integer(int32), parameter :: min_iter_count = 2
    integer(int32), parameter :: min_fun_count = 10
    integer(int32), parameter :: min_update_count = 1

    ! Local Variables
    logical :: stop
    integer(int32) :: m, n, actual, expected, flag
    real(real64), pointer :: w(:), pmax(:), pmin(:)
    real(real64), allocatable, target :: defaultWeights(:), maxparam(:), &
        minparam(:), JtWJ(:,:)
    type(iteration_controls) :: tol
    type(lm_solver_options) :: opt
    type(convergence_info) :: cInfo
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    type(convergence_info), target :: defaultinfo
    type(convergence_info), pointer :: inf
    
    ! Initialization
    stop = .false.
    m = size(x)
    n = size(params)
    if (present(info)) then
        inf => info
    else
        inf => defaultinfo
    end if
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (present(controls)) then
        tol = controls
    else
        call tol%set_to_default()
    end if
    if (present(settings)) then
        opt = settings
    else
        call opt%set_to_default()
    end if

    ! Input Checking
    if (size(y) /= m) then
        call report_array_size_error(errmgr, "nonlinear_least_squares", &
            "y", m, size(y))
        return
    end if
    if (size(ymod) /= m) then
        call report_array_size_error(errmgr, "nonlinear_least_squares", &
            "ymod", m, size(ymod))
        return
    end if
    if (size(resid) /= m) then
        call report_array_size_error(errmgr, "nonlinear_least_squares", &
            "resid", m, size(resid))
        return
    end if
    if (m < n) then
        call report_underdefined_error(errmgr, &
            "nonlinear_least_squares", n, m)
        return
    end if

    ! Tolerance Checking
    if (tol%gradient_tolerance < too_small) then
        call errmgr%report_error("nonlinear_least_squares", &
            "The gradient tolerance was found to be too small.", &
            FS_TOLERANCE_TOO_SMALL_ERROR)
        return
    end if
    if (tol%change_in_solution_tolerance < too_small) then
        call errmgr%report_error("nonlinear_least_squares", &
            "The change in solution tolerance was found to be too small.", &
            FS_TOLERANCE_TOO_SMALL_ERROR)
        return
    end if
    if (tol%residual_tolerance < too_small) then
        call errmgr%report_error("nonlinear_least_squares", &
            "The residual error tolerance was found to be too small.", &
            FS_TOLERANCE_TOO_SMALL_ERROR)
        return
    end if
    if (tol%iteration_improvement_tolerance < too_small) then
        call errmgr%report_error("nonlinear_least_squares", &
            "The iteration improvement tolerance was found to be too small.", &
            FS_TOLERANCE_TOO_SMALL_ERROR)
        return
    end if

    ! Iteration Count Checking
    if (tol%max_iteration_count < min_iter_count) then
        call report_iteration_count_error(errmgr, &
            "nonlinear_least_squares", &
            "Too few iterations were specified.", &
            min_iter_count)
        return
    end if
    if (tol%max_function_evaluations < min_fun_count) then
        call report_iteration_count_error(errmgr, &
            "nonlinear_least_squares", &
            "Too few function evaluations were specified.", &
            min_fun_count)
        return
    end if
    if (tol%max_iteration_between_updates < min_update_count) then
        call report_iteration_count_error(errmgr, &
            "nonlinear_least_squares", &
            "Too few iterations between updates were specified.", &
            min_update_count)
        return
    end if

    ! Optional Array Arguments (weights, parameter limits, etc.)
    if (present(weights)) then
        if (size(weights) < m) then
            call report_array_size_error(errmgr, &
                "nonlinear_least_squares", "weights", m, size(weights))
            return
        end if
        w(1:m) => weights(1:m)
    else
        allocate(defaultWeights(m), source = 1.0d0, stat = flag)
        if (flag /= 0) go to 50
        w(1:m) => defaultWeights(1:m)
    end if

    if (present(maxp)) then
        if (size(maxp) /= n) then
            call report_array_size_error(errmgr, &
                "nonlinear_least_squares", "maxp", n, size(maxp))
            return
        end if
        pmax(1:n) => maxp(1:n)
    else
        allocate(maxparam(n), source = huge(1.0d0), stat = flag)
        if (flag /= 0) go to 50
        pmax(1:n) => maxparam(1:n)
    end if

    if (present(minp)) then
        if (size(minp) /= n) then
            call report_array_size_error(errmgr, &
                "nonlinear_least_squares", "minp", n, size(minp))
            return
        end if
        pmin(1:n) => minp(1:n)
    else
        allocate(minparam(n), source = -huge(1.0d0), stat = flag)
        if (flag /= 0) go to 50
        pmin(1:n) => minparam(1:n)
    end if

    ! Local Memory Allocations
    allocate(JtWJ(n, n), stat = flag)
    if (flag /= 0) go to 50

    ! Process
    call lm_solve(fun, x, y, params, w, pmax, pmin, tol, opt, ymod, &
        resid, JtWJ, inf, stop, errmgr, status)

    ! Compute the covariance matrix
    if (present(stats) .or. present(cov)) then
        call mtx_inverse(JtWJ, err = errmgr)
        if (errmgr%has_error_occurred()) return
    end if

    ! Statistical Parameters
    if (present(stats)) then
        if (size(stats) /= n) then
            call report_array_size_error(errmgr, &
                "nonlinear_least_squares", "stats", n, size(stats))
            return
        end if

        ! Compute the statistics
        stats = calculate_regression_statistics(resid, params, JtWJ, &
            alpha, errmgr)
    end if

    ! Return the covariance matrix
    if (present(cov)) then
        if (size(cov, 1) /= n .or. size(cov, 2) /= n) then
            call report_matrix_size_error(errmgr, "nonlinear_least_squares", &
                "cov", n, n, size(cov, 1), size(cov, 2))
            return
        end if
        cov = JtWJ
    end if

    ! End
    return

    ! Memory Error Handler
50      continue
    call report_memory_error(errmgr, "nonlinear_least_squares", flag)
    return
end subroutine

! ******************************************************************************
! SETTINGS DEFAULTS
! ------------------------------------------------------------------------------
! Sets up default tolerances.
subroutine lm_set_default_tolerances(x)
    ! Arguments
    class(iteration_controls), intent(inout) :: x

    ! Set defaults
    x%max_iteration_count = 500
    x%max_function_evaluations = 5000
    x%max_iteration_between_updates = 10
    x%gradient_tolerance = 1.0d-8
    x%residual_tolerance = 0.5d-2
    x%change_in_solution_tolerance = 1.0d-6
    x%iteration_improvement_tolerance = 1.0d-1
end subroutine

! ------------------------------------------------------------------------------
! Sets up default solver settings.
subroutine lm_set_default_settings(x)
    ! Arguments
    class(lm_solver_options), intent(inout) :: x

    ! Set defaults
    x%method = FS_LEVENBERG_MARQUARDT_UPDATE
    x%finite_difference_step_size = sqrt(epsilon(1.0d0))
    x%damping_increase_factor = 11.0d0
    x%damping_decrease_factor = 9.0d0
end subroutine

! ******************************************************************************
! PRIVATE ROUTINES
! ------------------------------------------------------------------------------
! Computes the Jacobian matrix via a forward difference.
!
! Inputs:
! - fun: The function to evaluate
! - xdata: The independent coordinate data to fit (M-by-1)
! - params: The model parameters (N-by-1)
! - f0: The current model estimate (M-by-1)
! - step: The differentiation step size
!
! Outputs:
! - jac: The Jacobian matrix (M-by-N)
! - f1: A workspace array for the model output (M-by-1)
! - stop: A flag allowing the user to terminate model execution
! - work: A workspace array for the model parameters (N-by-1)
subroutine jacobian_finite_diff(fun, xdata, params, f0, jac, f1, &
    stop, step, work)
    ! Arguments
    procedure(regression_function), intent(in), pointer :: fun
    real(real64), intent(in) :: xdata(:), params(:)
    real(real64), intent(in) :: f0(:)
    real(real64), intent(out) :: jac(:,:)
    real(real64), intent(out) :: f1(:), work(:)
    logical, intent(out) :: stop
    real(real64), intent(in) :: step

    ! Local Variables
    integer(int32) :: i, n

    ! Initialization
    n = size(params)

    ! Cycle over each column of the Jacobian and calculate the derivative
    ! via a forward difference scheme
    !
    ! J(i,j) = df(i) / dx(j)
    work = params
    do i = 1, n
        work(i) = work(i) + step
        call fun(xdata, work, f1, stop)
        if (stop) return

        jac(:,i) = (f1 - f0) / step
        work(i) = params(i)
    end do
end subroutine

! ------------------------------------------------------------------------------
! Computes a rank-1 update to the Jacobian matrix
!
! Inputs:
! - pOld: previous set of parameters (N-by-1)
! - yOld: model evaluation at previous set of parameters (M-by-1)
! - jac: current Jacobian estimate (M-by-N)
! - p: current set of parameters (N-by-1)
! - y: model evaluation at current set of parameters (M-by-1)
! 
! Outputs:
! - jac: updated Jacobian matrix (M-by-N) (dy * dp**T + J)
! - dp: p - pOld (N-by-1)
! - dy: (y - yOld - J * dp) / (dp' * dp) (M-by-1)
subroutine broyden_update(pOld, yOld, jac, p, y, dp, dy)
    ! Arguments
    real(real64), intent(in) :: pOld(:), yOld(:), p(:), y(:)
    real(real64), intent(inout) :: jac(:,:)
    real(real64), intent(out) :: dp(:), dy(:)

    ! Local Variables
    real(real64) :: h2

    ! Process
    dp = p - pOld
    h2 = dot_product(dp, dp)
    dy = y - yOld - matmul(jac, dp)
    dy = dy / h2
    call rank1_update(1.0d0, dy, dp, jac)
end subroutine

! ------------------------------------------------------------------------------
! Updates the Levenberg-Marquardt matrix by either computing a new Jacobian
! matrix or performing a rank-1 update to the existing Jacobian matrix.
!
! Inputs:
! - fun: The function to evaluate
! - xdata: The independent coordinate data to fit (M-by-1)
! - ydata: The dependent coordinate data to fit (M-by-1)
! - pOld: previous set of parameters (N-by-1)
! - yOld: model evaluation at previous set of parameters (M-by-1)
! - dX2: The previous change in the Chi-squared criteria
! - jac: current Jacobian estimate (M-by-N)
! - p: current set of parameters (N-by-1)
! - weights: A weighting vector (M-by-1)
! - neval: Current number of function evaluations
! - update: Set to true to force an update of the Jacobian; else, set to
!       false to let the program choose based upon the change in the 
!       Chi-squared parameter.
! - step: The differentiation step size
!
! Outputs:
! - JtWJ: linearized Hessian matrix (inverse of the covariance matrix) (N-by-N)
! - JtWdy: linearized fitting vector (N-by-1)
! - X2: Updated Chi-squared criteria
! - yNew: model evaluated with parameters of p (M-by-1)
! - jac: updated Jacobian matrix (M-by-N)
! - neval: updated count of function evaluations
! - stop: A flag allowing the user to terminate model execution
! - work: A workspace array (N+M-by-1)
! - mwork: A workspace matrix (N-by-M)
! - update: Reset to false if a Jacobian evaluation was performed.
subroutine lm_matrix(fun, xdata, ydata, pOld, yOld, dX2, jac, p, weights, &
    neval, update, step, JtWJ, JtWdy, X2, yNew, stop, work, mwork)
    ! Arguments
    procedure(regression_function), pointer :: fun
    real(real64), intent(in) :: xdata(:), ydata(:), pOld(:), yOld(:), &
        p(:), weights(:)
    real(real64), intent(in) :: dX2, step
    real(real64), intent(inout) :: jac(:,:)
    integer(int32), intent(inout) :: neval
    logical, intent(inout) :: update
    real(real64), intent(out) :: JtWJ(:,:), JtWdy(:)
    real(real64), intent(out) :: X2, mwork(:,:), yNew(:)
    logical, intent(out) :: stop
    real(real64), intent(out), target :: work(:)

    ! Local Variables
    integer(int32) :: m, n
    real(real64), pointer :: w1(:), w2(:)

    ! Initialization
    m = size(xdata)
    n = size(p)
    w1(1:m) => work(1:m)
    w2(1:n) => work(m+1:n+m)

    ! Perform the next function evaluation
    call fun(xdata, p, yNew, stop)
    neval = neval + 1
    if (stop) return

    ! Update or recompute the Jacobian matrix
    if (dX2 > 0 .or. update) then
        ! Recompute the Jacobian
        call jacobian_finite_diff(fun, xdata, p, yNew, jac, w1, &
            stop, step, w2)
        neval = neval + n
        if (stop) return
        update = .false.
    else
        ! Simply perform a rank-1 update to the Jacobian
        call broyden_update(pOld, yOld, jac, p, yNew, w2, w1)
    end if

    ! Update the Chi-squared estimate
    w1 = ydata - yNew
    X2 = dot_product(w1, w1 * weights)

    ! Compute J**T * (W .* dY)
    w1 = w1 * weights
    call mtx_mult(.true., 1.0d0, jac, w1, 0.0d0, JtWdy)

    ! Update the Hessian
    ! First: J**T * W = MWORK
    ! Second: (J**T * W) * J
    call diag_mtx_mult(.false., .true., 1.0d0, weights, jac, 0.0d0, mwork)
    call mtx_mult(.false., .false., 1.0d0, mwork, jac, 0.0d0, JtWJ)
end subroutine

! ------------------------------------------------------------------------------
! Performs a single iteration of the Levenberg-Marquardt algorithm.
!
! Inputs:
! - fun: The function to evaluate
! - xdata: The independent coordinate data to fit (M-by-1)
! - ydata: The dependent coordinate data to fit (M-by-1)
! - p: current set of parameters (N-by-1)
! - neval: current number of function evaluations
! - niter: current iteration number
! - update: set to 1 to use Marquardt's modification; else, 
! - step: the differentiation step size
! - lambda: LM damping parameter
! - maxP: maximum limits on the parameters.  Use huge() or larger for no constraints (N-by-1)
! - minP: minimum limits on the parameters.  Use -huge() or smaller for no constraints (N-by-1)
! - weights: a weighting vector (M-by-1)
! - JtWJ: linearized Hessian matrix (inverse of the covariance matrix) (N-by-N)
! - JtWdy: linearized fitting vector (N-by-1)
!
! Outputs:
! - JtWJ: overwritten LU factorization of the original matrix (N-by-N)
! - h: The new estimate of the change in parameter (N-by-1)
! - pNew: The new parameter estimates (N-by-1)
! - deltaY: The new difference between data and model (M-by-1)
! - yNew: model evaluated with parameters of pNew (M-by-1)
! - neval: updated count of function evaluations
! - niter: updated current iteration number
! - X2: updated Chi-squared criteria
! - stop: A flag allowing the user to terminate model execution
! - iwork: A workspace array (N-by-1)
! - err: An error handling mechanism
subroutine lm_iter(fun, xdata, ydata, p, neval, niter, update, lambda, &
    maxP, minP, weights, JtWJ, JtWdy, h, pNew, deltaY, yNew, X2, X2Old, &
    alpha, stop, iwork, err, status)
    ! Arguments
    procedure(regression_function), pointer :: fun
    real(real64), intent(in) :: xdata(:), ydata(:), p(:), maxP(:), &
        minP(:), weights(:), JtWdy(:)
    real(real64), intent(in) :: lambda, X2Old
    integer(int32), intent(inout) :: neval, niter
    integer(int32), intent(in) :: update
    real(real64), intent(inout) :: JtWJ(:,:)
    real(real64), intent(out) :: h(:), pNew(:), deltaY(:), yNew(:)
    real(real64), intent(out) :: X2, alpha
    logical, intent(out) :: stop
    integer(int32), intent(out) :: iwork(:)
    class(errors), intent(inout) :: err
    procedure(iteration_update), intent(in), pointer, optional :: status

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: dpJh

    ! Initialization
    n = size(p)

    ! Increment the iteration counter
    niter = niter + 1

    ! Solve the linear system to determine the change in parameters
    ! A is N-by-N and is stored in JtWJ
    ! b is N-by-1
    if (update == FS_LEVENBERG_MARQUARDT_UPDATE) then
        ! Compute: h = A \ b
        ! A = J**T * W * J + lambda * diag(J**T * W * J)
        ! b = J**T * W * dy
        do i = 1, n
            JtWJ(i,i) = JtWJ(i,i) * (1.0d0 + lambda)
            h(i) = JtWdy(i)
        end do
    else
        ! Compute: h = A \ b
        ! A = J**T * W * J + lambda * I
        ! b = J**T * W * dy
        do i = 1, n
            JtWJ(i,i) = JtWJ(i,i) + lambda
            h(i) = JtWdy(i)
        end do
    end if
    call lu_factor(JtWJ, iwork, err)        ! overwrites JtWJ with [L\U]
    if (err%has_error_occurred()) return    ! if JtWJ is singular
    call solve_lu(JtWJ, iwork, h)           ! solution stored in h

    ! Compute the new attempted solution, and apply any constraints
    do i = 1, n
        pNew(i) = min(max(minP(i), h(i) + p(i)), maxP(i))
    end do

    ! Update the residual error
    call fun(xdata, pNew, yNew, stop)
    neval = neval + 1
    deltaY = ydata - yNew
    if (stop) return

    ! Update the Chi-squared estimate
    X2 = dot_product(deltaY, deltaY * weights)

    ! Perform a quadratic line update in the H direction, if necessary
    if (update == FS_QUADRATIC_UPDATE) then
        dpJh = dot_product(JtWdy, h)
        alpha = abs(dpJh / (0.5d0 * (X2 - X2Old) + 2.0d0 * dpJh))
        h = alpha * h

        do i = 1, n
            pNew(i) = min(max(minP(i), p(i) + h(i)), maxP(i))
        end do

        call fun(xdata, pNew, yNew, stop)
        if (stop) return
        neval = neval + 1
        deltaY = ydata - yNew
        X2 = dot_product(deltaY, deltaY * weights)
    end if

    ! Update the status of the iteration, if needed
    if (present(status)) then
        call status(niter, yNew, deltaY, pNew, h)
    end if
end subroutine

! ------------------------------------------------------------------------------
! A Levenberg-Marquardt solver.
!
! Inputs:
! - fun: The function to evaluate
! - xdata: The independent coordinate data to fit (M-by-1)
! - ydata: The dependent coordinate data to fit (M-by-1)
! - p: current set of parameters (N-by-1)
! - weights: a weighting vector (M-by-1)
! - maxP: maximum limits on the parameters.  Use huge() or larger for no constraints (N-by-1)
! - minP: minimum limits on the parameters.  Use -huge() or smaller for no constraints (N-by-1)
! - controls: an iteration_controls instance containing solution tolerances
!
! Outputs:
! - p: solution (N-by-1)
! - y: model results at p (M-by-1)
! - resid: residual (ydata - y) (M-by-1)
! - JtWJ: linearized Hessian matrix (inverse of the covariance matrix) (N-by-N)
! - opt: a convergence_info object containing information regarding 
!       convergence of the iteration
! - stop: A flag allowing the user to terminate model execution
! - err: An error handling object
subroutine lm_solve(fun, xdata, ydata, p, weights, maxP, minP, controls, &
    opt, y, resid, JtWJ, info, stop, err, status)
    ! Arguments
    procedure(regression_function), intent(in), pointer :: fun
    real(real64), intent(in) :: xdata(:), ydata(:), weights(:), maxP(:), &
        minP(:)
    real(real64), intent(inout) :: p(:)
    class(iteration_controls), intent(in) :: controls
    class(lm_solver_options), intent(in) :: opt
    real(real64), intent(out) :: y(:), resid(:), JtWJ(:,:)
    class(convergence_info), intent(out) :: info
    logical, intent(out) :: stop
    class(errors), intent(inout) :: err
    procedure(iteration_update), intent(in), pointer, optional :: status

    ! Local Variables
    logical :: update
    integer(int32) :: i, m, n, dof, flag, neval, niter, nupdate
    real(real64) :: dX2, X2, X2Old, X2Try, lambda, alpha, nu, step
    real(real64), allocatable :: pOld(:), yOld(:), J(:,:), JtWdy(:), &
        work(:), mwork(:,:), pTry(:), yTemp(:), JtWJc(:,:), h(:)
    integer(int32), allocatable :: iwork(:)
    character(len = :), allocatable :: errmsg

    ! Initialization
    update = .true.
    m = size(xdata)
    n = size(p)
    dof = m - n
    niter = 0
    step = opt%finite_difference_step_size
    stop = .false.
    info%user_requested_stop = .false.
    nupdate = 0

    ! Local Memory Allocation
    allocate(pOld(n), source = 0.0d0, stat = flag)
    if (flag == 0) allocate(yOld(m), source = 0.0d0, stat = flag)
    if (flag == 0) allocate(J(m, n), stat = flag)
    if (flag == 0) allocate(JtWdy(n), stat = flag)
    if (flag == 0) allocate(work(m + n), stat = flag)
    if (flag == 0) allocate(mwork(n, m), stat = flag)
    if (flag == 0) allocate(pTry(n), stat = flag)
    if (flag == 0) allocate(h(n), stat = flag)
    if (flag == 0) allocate(yTemp(m), stat = flag)
    if (flag == 0) allocate(JtWJc(n, n), stat = flag)
    if (flag == 0) allocate(iwork(n), stat = flag)
    if (flag /= 0) go to 10

    ! Perform an initial function evaluation
    call fun(xdata, p, y, stop)
    neval = 1

    ! Evaluate the problem matrices
    call lm_matrix(fun, xdata, ydata, pOld, yOld, 1.0d0, J, p, weights, &
        neval, update, step, JtWJ, JtWdy, X2, y, stop, work, mwork)
    if (stop) go to 5
    X2Old = X2
    JtWJc = JtWJ

    ! Determine an initial value for lambda
    if (opt%method == FS_LEVENBERG_MARQUARDT_UPDATE) then
        lambda = 1.0d-2
    else
        call extract_diagonal(JtWJ, work(1:n))
        lambda = 1.0d-2 * maxval(work(1:n))
        nu = 2.0d0
    end if

    ! Main Loop
    main : do while (niter < controls%max_iteration_count)
        ! Compute the linear solution at the current solution estimate and
        ! update the new parameter estimates
        call lm_iter(fun, xdata, ydata, p, neval, niter, opt%method, &
            lambda, maxP, minP, weights, JtWJc, JtWdy, h, pTry, resid, &
            yTemp, X2Try, X2Old, alpha, stop, iwork, err, status)
        if (stop) go to 5
        if (err%has_error_occurred()) return

        ! Update the Chi-squared estimate, update the damping parameter
        ! lambda, and, if necessary, update the matrices
        call lm_update(fun, xdata, ydata, pOld, p, pTry, yOld, y, h, dX2, &
            X2Old, X2, X2Try, lambda, alpha, nu, JtWdy, JtWJ, J, weights, &
            niter, neval, update, step, work, mwork, controls, opt, stop)
        if (stop) go to 5
        JtWJc = JtWJ

        ! Determine the matrix update scheme
        nupdate = nupdate + 1
        if (opt%method == FS_QUADRATIC_UPDATE) then
            update = mod(niter, 2 * n) > 0
        else if (nupdate >= controls%max_iteration_between_updates) then
            update = .true.
            nupdate = 0
        end if

        ! Test for convergence
        if (lm_check_convergence(controls, dof, resid, niter, neval, &
            JtWdy, h, p, X2, info)) &
        then
            exit main
        end if
    end do main

    ! End
    return

    ! User Requested End
5       continue
    info%user_requested_stop = .true.
    return

    ! Memory Error Handling
10      continue
    allocate(character(len = 512) :: errmsg)
    write(errmsg, 100) "Memory allocation error code ", flag, "."
    call err%report_error("lm_solve", &
        trim(errmsg), FS_MEMORY_ERROR)
    return

    ! Formatting
100     format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
!
subroutine lm_update(fun, xdata, ydata, pOld, p, pTry, yOld, y, h, dX2, &
    X2old, X2, X2try, lambda, alpha, nu, JtWdy, JtWJ, J, weights, niter, &
    neval, update, step, work, mwork, controls, opt, stop)
    ! Arguments
    procedure(regression_function), intent(in), pointer :: fun
    real(real64), intent(in) :: xdata(:), ydata(:), X2try, h(:), step, &
        pTry(:), weights(:), alpha
    real(real64), intent(inout) :: pOld(:), p(:), yOld(:), y(:), lambda, &
        JtWdy(:), dX2, X2, X2old, JtWJ(:,:), J(:,:), nu
    real(real64), intent(out) :: work(:), mwork(:,:)
    integer(int32), intent(in) :: niter
    integer(int32), intent(inout) :: neval
    logical, intent(inout) :: update
    class(iteration_controls), intent(in) :: controls
    class(lm_solver_options), intent(in) :: opt
    logical, intent(out) :: stop

    ! Local Variables
    integer(int32) :: n
    real(real64) :: rho

    ! Initialization
    n = size(p)

    ! Process
    if (opt%method == FS_LEVENBERG_MARQUARDT_UPDATE) then
        call extract_diagonal(JtWJ, work(1:n))
        work(1:n) = lambda * work(1:n) * h + JtWdy
    else
        work(1:n) = lambda * h + JtWdy
    end if
    rho = (X2 - X2try) / abs(dot_product(h, work(1:n)))
    if (rho > controls%iteration_improvement_tolerance) then
        ! Things are getting better at an acceptable rate
        dX2 = X2 - X2old
        X2old = X2
        pOld = p
        yOld = y
        p = pTry

        ! Recompute the matrices
        call lm_matrix(fun, xdata, ydata, pOld, yOld, dX2, J, p, weights, &
            neval, update, step, JtWJ, JtWdy, X2, y, stop, work, mwork)
        if (stop) return

        ! Decrease lambda
        select case (opt%method)
        case (FS_LEVENBERG_MARQUARDT_UPDATE)
            lambda = max(lambda / opt%damping_decrease_factor, 1.0d-7)
        case (FS_QUADRATIC_UPDATE)
            lambda = max(lambda / (1.0d0 + alpha), 1.0d-7)
        case (FS_NIELSEN_UPDATE)
            lambda = lambda * max(1.0d0 / 3.0d0, &
                1.0d0 - (2.0d0 * rho - 1.0d0**3))
            nu = 2.0d0
        end select
    else
        ! The iteration is not improving in a satisfactory manner
        X2 = X2old
        if (mod(niter, 2 * n) /= 0) then
            call lm_matrix(fun, xdata, ydata, pOld, yOld, -1.0d0, J, p, &
                weights, neval, update, step, JtWJ, JtWdy, dX2, y, stop, &
                work, mwork)
            if (stop) return
        end if

        ! Increase lambda
        select case (opt%method)
        case (FS_LEVENBERG_MARQUARDT_UPDATE)
            lambda = min(lambda * opt%damping_increase_factor, 1.0d7)
        case (FS_QUADRATIC_UPDATE)
            lambda = lambda + abs((X2try - X2) / 2.0d0 / alpha)
        case (FS_NIELSEN_UPDATE)
            lambda = lambda * nu
            nu = 2.0d0 * nu
        end select
    end if
end subroutine

! ------------------------------------------------------------------------------
! Checks the Levenberg-Marquardt solution against the convergence criteria.
!
! Inputs:
! - controls: the solution controls and convergence criteria
! - dof: the statistical degrees of freedom of the system (M - N)
! - resid: the residual error (M-by-1)
! - niter: the number of iterations
! - neval: the number of function evaluations
! - JtWdy: linearized fitting vector (N-by-1)
! - h: the change in parameter (solution) values (N-by-1)
! - p: the parameter (solution) values (N-by-1)
! - X2: the Chi-squared estimate
!
! Outputs:
! - info: The convergence information.
! - rst: True if convergence was achieved; else, false.
function lm_check_convergence(controls, dof, resid, niter, neval, &
    JtWdy, h, p, X2, info) result(rst)
    ! Arguments
    class(iteration_controls), intent(in) :: controls
    real(real64), intent(in) :: resid(:), JtWdy(:), h(:), p(:), X2
    integer(int32), intent(in) :: dof, niter, neval
    class(convergence_info), intent(out) :: info
    logical :: rst

    ! Initialization
    rst = .false.

    ! Iteration Checks
    info%iteration_count = niter
    if (niter >= controls%max_iteration_count) then
        info%reach_iteration_limit = .true.
        rst = .true.
    else
        info%reach_iteration_limit = .false.
    end if

    info%function_evaluation_count = neval
    if (neval >= controls%max_function_evaluations) then
        info%reach_function_evaluation_limit = .true.
        rst = .true.
    else
        info%reach_function_evaluation_limit = .false.
    end if

    info%gradient_value = maxval(abs(JtWdy))
    if (info%gradient_value < controls%gradient_tolerance .and. niter > 2) &
    then
        info%converge_on_gradient = .true.
        rst = .true.
    else
        info%converge_on_gradient = .false.
    end if

    info%solution_change_value = maxval(abs(h) / (abs(p) + 1.0d-12))
    if (info%solution_change_value < &
        controls%change_in_solution_tolerance .and. niter > 2) &
    then
        info%converge_on_solution_change = .true.
        rst = .true.
    else
        info%converge_on_solution_change = .false.
    end if

    info%residual_value = X2 / dof
    if (info%residual_value < controls%residual_tolerance .and. niter > 2) &
    then
        info%converge_on_residual_parameter = .true.
        rst = .true.
    else
        info%converge_on_residual_parameter = .false.
    end if
end function

! ------------------------------------------------------------------------------
end module
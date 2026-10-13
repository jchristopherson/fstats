module fstats_levenberg_marquardt
    use iso_fortran_env
    use ieee_arithmetic
    use fstats_errors
    use fstats_regression
    use linalg
    implicit none
    private
    public :: lm_solver_options
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

contains
! ------------------------------------------------------------------------------
subroutine nonlinear_least_squares(fun, x, y, params, ymod, &
    resid, weights, maxp, minp, stats, alpha, controls, settings, info, &
    status, cov, args)
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
        !! the covariance matrix, computed from an SVD of the final weighted
        !! Jacobian. Rank diagnostics are available in info when requested.
    class(*), intent(inout), optional :: args
        !! An optional argument allowing the passing in/out of data for the
        !! [[fun]] routine.

    ! Parameters
    real(real64), parameter :: too_small = 1.0d-14
    integer(int32), parameter :: min_iter_count = 2
    integer(int32), parameter :: min_fun_count = 10
    integer(int32), parameter :: min_update_count = 1

    ! Local Variables
    logical :: stop
    integer(int32) :: m, n, actual, expected
    real(real64), pointer :: w(:), pmax(:), pmin(:)
    real(real64), allocatable, target :: defaultWeights(:), maxparam(:), &
        minparam(:), JtWJ(:,:)
    real(real64), allocatable, dimension(:,:) :: final_jacobian, final_covariance
    type(iteration_controls) :: tol
    type(lm_solver_options) :: opt
    type(convergence_info) :: cInfo
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
    if (size(y) /= m) error stop FS_ARRAY_SIZE_ERROR
    if (size(ymod) /= m) error stop FS_ARRAY_SIZE_ERROR
    if (size(resid) /= m) error stop FS_ARRAY_SIZE_ERROR
    if (m < n) error stop FS_UNDERDEFINED_PROBLEM_ERROR

    ! Tolerance Checking
    if (tol%gradient_tolerance < too_small) error stop FS_TOLERANCE_TOO_SMALL_ERROR
    if (tol%change_in_solution_tolerance < too_small) error stop FS_TOLERANCE_TOO_SMALL_ERROR
    if (tol%residual_tolerance < too_small) error stop FS_TOLERANCE_TOO_SMALL_ERROR
    if (tol%iteration_improvement_tolerance < too_small) error stop FS_TOLERANCE_TOO_SMALL_ERROR

    ! Iteration Count Checking
    if (tol%max_iteration_count < min_iter_count) error stop FS_TOO_FEW_ITERATION_ERROR
    if (tol%max_function_evaluations < min_fun_count) error stop FS_TOO_FEW_ITERATION_ERROR
    if (tol%max_iteration_between_updates < min_update_count) error stop FS_TOO_FEW_ITERATION_ERROR

    ! Optional Array Arguments (weights, parameter limits, etc.)
    if (present(weights)) then
        if (size(weights) < m) error stop FS_ARRAY_SIZE_ERROR
        w(1:m) => weights(1:m)
    else
        allocate(defaultWeights(m), source = 1.0d0)
        w(1:m) => defaultWeights(1:m)
    end if
    if (.not.all(ieee_is_finite(w)) .or. any(w < 0.0d0)) error stop FS_INVALID_INPUT_ERROR

    if (present(maxp)) then
        if (size(maxp) /= n) error stop FS_ARRAY_SIZE_ERROR
        pmax(1:n) => maxp(1:n)
    else
        allocate(maxparam(n), source = huge(1.0d0))
        pmax(1:n) => maxparam(1:n)
    end if

    if (present(minp)) then
        if (size(minp) /= n) error stop FS_ARRAY_SIZE_ERROR
        pmin(1:n) => minp(1:n)
    else
        allocate(minparam(n), source = -huge(1.0d0))
        pmin(1:n) => minparam(1:n)
    end if

    ! Local Memory Allocations
    allocate(JtWJ(n, n))

    ! Process
    call lm_solve(fun, x, y, params, w, pmax, pmin, tol, opt, ymod, &
        resid, JtWJ, inf, stop, status, args = args)

    ! Compute the covariance matrix
    if (present(stats) .or. present(cov)) then
        allocate(final_jacobian(m, n))
        call jacobian(fun, x, params, final_jacobian, stop, f0 = ymod, &
            step = opt%finite_difference_step_size, args = args)
        inf%function_evaluation_count = inf%function_evaluation_count + n
        if (stop) then
            inf%user_requested_stop = .true.
            JtWJ = ieee_value(0.0d0, ieee_quiet_nan)
        else
            final_jacobian = final_jacobian * spread(sqrt(w), 2, n)
            call regression_covariance(final_jacobian, final_covariance, &
                inf%covariance_rank, inf%covariance_condition_number)
            JtWJ = final_covariance
        end if
    end if

    ! Statistical Parameters
    if (present(stats)) then
        if (size(stats) /= n) error stop FS_ARRAY_SIZE_ERROR

        ! Compute the statistics
        stats = calculate_regression_statistics(resid, params, JtWJ, alpha)
    end if

    ! Return the covariance matrix
    if (present(cov)) then
        if (size(cov, 1) /= n .or. size(cov, 2) /= n) error stop FS_MATRIX_SIZE_ERROR
        cov = JtWJ
    end if
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
    neval, update, step, JtWJ, JtWdy, X2, yNew, stop, work, mwork, args)
    ! Arguments
    procedure(regression_function), pointer :: fun
    real(real64), intent(in), dimension(:) :: xdata, ydata, pOld, yOld, &
        p, weights
    real(real64), intent(in) :: dX2, step
    real(real64), intent(inout), dimension(:,:) :: jac
    integer(int32), intent(inout) :: neval
    logical, intent(inout) :: update
    real(real64), intent(out), dimension(:,:) :: JtWJ
    real(real64), intent(out), dimension(:) :: JtWdy
    real(real64), intent(out) :: X2
    real(real64), intent(out), dimension(:,:) :: mwork
    real(real64), intent(out), dimension(:) :: yNew
    logical, intent(out) :: stop
    real(real64), intent(out), target, dimension(:) :: work
    class(*), intent(inout), optional :: args

    ! Local Variables
    integer(int32) :: m, n
    real(real64), pointer, dimension(:) :: w1, w2

    ! Initialization
    m = size(xdata)
    n = size(p)
    w1(1:m) => work(1:m)
    w2(1:n) => work(m+1:n+m)

    ! Perform the next function evaluation
    call fun(xdata, p, yNew, stop, args = args)
    neval = neval + 1
    if (stop) return

    ! Update or recompute the Jacobian matrix
    if (dX2 > 0 .or. update) then
        ! Recompute the Jacobian
        call jacobian(fun, xdata, p, jac, stop, f0 = yNew, f1 = w1, &
            step = step, args = args)
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
! - h: The new estimate of the change in parameter (N-by-1)
! - pNew: The new parameter estimates (N-by-1)
! - deltaY: The new difference between data and model (M-by-1)
! - yNew: model evaluated with parameters of pNew (M-by-1)
! - neval: updated count of function evaluations
! - niter: updated current iteration number
! - X2: updated Chi-squared criteria
! - stop: A flag allowing the user to terminate model execution
subroutine lm_iter(fun, xdata, ydata, p, neval, niter, update, lambda, &
    maxP, minP, weights, JtWJ, JtWdy, h, pNew, deltaY, yNew, X2, X2Old, &
    alpha, stop, status, args)
    ! Arguments
    procedure(regression_function), pointer :: fun
    real(real64), intent(in) :: xdata(:), ydata(:), p(:), maxP(:), &
        minP(:), weights(:), JtWdy(:)
    real(real64), intent(in) :: lambda, X2Old
    integer(int32), intent(inout) :: neval, niter
    integer(int32), intent(in) :: update
    real(real64), intent(in) :: JtWJ(:,:)
    real(real64), intent(out) :: h(:), pNew(:), deltaY(:), yNew(:)
    real(real64), intent(out) :: X2, alpha
    logical, intent(out) :: stop
    procedure(iteration_update), intent(in), pointer, optional :: status
    class(*), intent(inout), optional :: args

    ! Local Variables
    integer(int32) :: i, n
    integer(int32), allocatable, dimension(:) :: iwork
    real(real64) :: dpJh
    real(real64), allocatable, dimension(:,:) :: a, lu

    ! Initialization
    n = size(p)
    allocate(a(n,n), source = JtWJ)

    ! Increment the iteration counter
    niter = niter + 1

    ! Solve the linear system to determine the change in parameters
    ! b is N-by-1
    if (update == FS_LEVENBERG_MARQUARDT_UPDATE) then
        ! Compute: h = A \ b
        ! A = J**T * W * J + lambda * diag(J**T * W * J)
        ! b = J**T * W * dy
        do i = 1, n
            a(i,i) = a(i,i) * (1.0d0 + lambda)
        end do
    else
        ! Compute: h = A \ b
        ! A = J**T * W * J + lambda * I
        ! b = J**T * W * dy
        do i = 1, n
            a(i,i) = a(i,i) + lambda
        end do
    end if
    call lu_factor(a, ipvt = iwork, lu = lu)
    h = solve_lu(lu, iwork, JtWdy)

    ! Compute the new attempted solution, and apply any constraints
    do i = 1, n
        pNew(i) = min(max(minP(i), h(i) + p(i)), maxP(i))
    end do

    ! Update the residual error
    call fun(xdata, pNew, yNew, stop, args = args)
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

        call fun(xdata, pNew, yNew, stop, args = args)
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
subroutine lm_solve(fun, xdata, ydata, p, weights, maxP, minP, controls, &
    opt, y, resid, JtWJ, info, stop, status, args)
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
    procedure(iteration_update), intent(in), pointer, optional :: status
    class(*), intent(inout), optional :: args

    ! Local Variables
    logical :: update
    integer(int32) :: i, m, n, dof, flag, neval, niter, nupdate
    real(real64) :: dX2, X2, X2Old, X2Try, lambda, alpha, nu, step
    real(real64), allocatable :: pOld(:), yOld(:), J(:,:), JtWdy(:), &
        work(:), mwork(:,:), pTry(:), yTemp(:), JtWJc(:,:), h(:)

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
    allocate(pOld(n), source = 0.0d0)
    allocate(yOld(m), source = 0.0d0)
    allocate( &
        J(m, n), &
        JtWdy(n), &
        work(m + n), &
        mwork(n, m), &
        pTry(n), &
        h(n), &
        yTemp(m), &
        JtWJc(n, n) &
    )

    ! Perform an initial function evaluation
    call fun(xdata, p, y, stop, args = args)
    neval = 1

    ! Evaluate the problem matrices
    call lm_matrix(fun, xdata, ydata, pOld, yOld, 1.0d0, J, p, weights, &
        neval, update, step, JtWJ, JtWdy, X2, y, stop, work, mwork, args = args)
    if (stop) go to 5
    X2Old = X2
    JtWJc = JtWJ

    ! Determine an initial value for lambda
    if (opt%method == FS_LEVENBERG_MARQUARDT_UPDATE) then
        lambda = 1.0d-2
    else
        work(1:n) = extract_diagonal(JtWJ)
        lambda = 1.0d-2 * maxval(work(1:n))
        nu = 2.0d0
    end if

    ! Main Loop
    main : do while (niter < controls%max_iteration_count)
        ! Compute the linear solution at the current solution estimate and
        ! update the new parameter estimates
        call lm_iter(fun, xdata, ydata, p, neval, niter, opt%method, &
            lambda, maxP, minP, weights, JtWJc, JtWdy, h, pTry, resid, &
            yTemp, X2Try, X2Old, alpha, stop, status, args = args)
        if (stop) go to 5

        ! Update the Chi-squared estimate, update the damping parameter
        ! lambda, and, if necessary, update the matrices
        call lm_update(fun, xdata, ydata, pOld, p, pTry, yOld, y, h, dX2, &
            X2Old, X2, X2Try, lambda, alpha, nu, JtWdy, JtWJ, J, weights, &
            niter, neval, update, step, work, mwork, controls, opt, stop, &
            args = args)
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
end subroutine

! ------------------------------------------------------------------------------
!
subroutine lm_update(fun, xdata, ydata, pOld, p, pTry, yOld, y, h, dX2, &
    X2old, X2, X2try, lambda, alpha, nu, JtWdy, JtWJ, J, weights, niter, &
    neval, update, step, work, mwork, controls, opt, stop, args)
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
    class(*), intent(inout), optional :: args

    ! Local Variables
    integer(int32) :: n
    real(real64) :: rho

    ! Initialization
    n = size(p)

    ! Process
    if (opt%method == FS_LEVENBERG_MARQUARDT_UPDATE) then
        work(1:n) = extract_diagonal(JtWJ)
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
            neval, update, step, JtWJ, JtWdy, X2, y, stop, work, mwork, &
            args = args)
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
                work, mwork, args = args)
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

! ******************************************************************************
! SETTINGS DEFAULTS
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

! ------------------------------------------------------------------------------
end module
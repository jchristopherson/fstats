submodule (fstats) levenberg_marquardt
! REFERENCES:
! 1. https://people.duke.edu/~hpgavin/ExperimentalSystems/lm.pdf
    use linalg, only : mtx_inverse, lu_factor, mtx_mult, solve_lu, &
        diag_mtx_mult, rank1_update
    use fstats_errors
    implicit none
contains
! ------------------------------------------------------------------------------
    module subroutine regression_jacobian_1(fun, xdata, params, &
        jac, stop, f0, f1, step, err)
        ! Arguments
        procedure(regression_function), intent(in), pointer :: fun
        real(real64), intent(in) :: xdata(:), params(:)
        real(real64), intent(out) :: jac(:,:)
        logical, intent(out) :: stop
        real(real64), intent(in), optional, target :: f0(:)
        real(real64), intent(out), optional, target :: f1(:)
        real(real64), intent(in), optional :: step
        class(errors), intent(inout), optional, target :: err

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
            call report_matrix_size_error(errmgr, "regression_jacobian_1", &
                "JAC", m, n, size(jac, 1), size(jac, 2))
            return
        end if
        if (present(f0)) then
            ! Check Size
            if (size(f0) /= m) then
                call report_array_size_error(errmgr, "regression_jacobian_1", &
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
                call report_array_size_error(errmgr, "regression_jacobian_1", &
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
20      continue
        call report_memory_error(errmgr, "regression_jacobian_1", flag)
        return
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine nonlinear_least_squares_1(fun, x, y, params, ymod, &
        resid, weights, maxp, minp, stats, alpha, controls, settings, info, &
        status, err)
        ! Arguments
        procedure(regression_function), intent(in), pointer :: fun
        real(real64), intent(in) :: x(:), y(:)
        real(real64), intent(inout) :: params(:)
        real(real64), intent(out) :: ymod(:), resid(:)
        real(real64), intent(in), optional, target :: weights(:), maxp(:), &
            minp(:)
        type(regression_statistics), intent(out), optional :: stats(:)
        real(real64), intent(in), optional :: alpha
        type(iteration_controls), intent(in), optional :: controls
        type(lm_solver_options), intent(in), optional :: settings
        type(convergence_info), intent(out), optional, target :: info
        procedure(iteration_update), intent(in), pointer, optional :: status
        class(errors), intent(inout), optional, target :: err

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
            call report_array_size_error(errmgr, "nonlinear_least_squares_1", &
                "y", m, size(y))
            return
        end if
        if (size(ymod) /= m) then
            call report_array_size_error(errmgr, "nonlinear_least_squares_1", &
                "ymod", m, size(ymod))
            return
        end if
        if (size(resid) /= m) then
            call report_array_size_error(errmgr, "nonlinear_least_squares_1", &
                "resid", m, size(resid))
            return
        end if
        if (m < n) then
            call report_underdefined_error(errmgr, &
                "nonlinear_least_squares_1", n, m)
            return
        end if

        ! Tolerance Checking
        if (tol%gradient_tolerance < too_small) then
            call errmgr%report_error("nonlinear_least_squares_1", &
                "The gradient tolerance was found to be too small.", &
                FS_TOLERANCE_TOO_SMALL_ERROR)
            return
        end if
        if (tol%change_in_solution_tolerance < too_small) then
            call errmgr%report_error("nonlinear_least_squares_1", &
                "The change in solution tolerance was found to be too small.", &
                FS_TOLERANCE_TOO_SMALL_ERROR)
            return
        end if
        if (tol%residual_tolerance < too_small) then
            call errmgr%report_error("nonlinear_least_squares_1", &
                "The residual error tolerance was found to be too small.", &
                FS_TOLERANCE_TOO_SMALL_ERROR)
            return
        end if
        if (tol%iteration_improvement_tolerance < too_small) then
            call errmgr%report_error("nonlinear_least_squares_1", &
                "The iteration improvement tolerance was found to be too small.", &
                FS_TOLERANCE_TOO_SMALL_ERROR)
            return
        end if

        ! Iteration Count Checking
        if (tol%max_iteration_count < min_iter_count) then
            call report_iteration_count_error(errmgr, &
                "nonlinear_least_squares_1", &
                "Too few iterations were specified.", &
                min_iter_count)
            return
        end if
        if (tol%max_function_evaluations < min_fun_count) then
            call report_iteration_count_error(errmgr, &
                "nonlinear_least_squares_1", &
                "Too few function evaluations were specified.", &
                min_fun_count)
            return
        end if
        if (tol%max_iteration_between_updates < min_update_count) then
            call report_iteration_count_error(errmgr, &
                "nonlinear_least_squares_1", &
                "Too few iterations between updates were specified.", &
                min_update_count)
            return
        end if

        ! Optional Array Arguments (weights, parameter limits, etc.)
        if (present(weights)) then
            if (size(weights) < m) then
                call report_array_size_error(errmgr, &
                    "nonlinear_least_squares_1", "weights", m, size(weights))
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
                    "nonlinear_least_squares_1", "maxp", n, size(maxp))
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
                    "nonlinear_least_squares_1", "minp", n, size(minp))
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

        ! Statistical Parameters
        if (present(stats)) then
            if (size(stats) /= n) then
                call report_array_size_error(errmgr, &
                    "nonlinear_least_squares_1", "stats", n, size(stats))
                return
            end if

            ! Compute the covariance matrix
            call mtx_inverse(JtWJ, err = errmgr)
            if (errmgr%has_error_occurred()) return

            ! Compute the statistics
            stats = calculate_regression_statistics(resid, params, JtWJ, &
                alpha, errmgr)
        end if

        ! End
        return

        ! Memory Error Handler
50      continue
        call report_memory_error(errmgr, "nonlinear_least_squares_1", flag)
        return
    end subroutine

! ******************************************************************************
! SETTINGS DEFAULTS
! ------------------------------------------------------------------------------
    ! Sets up default tolerances.
    module subroutine lm_set_default_tolerances(x)
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
    module subroutine lm_set_default_settings(x)
        ! Arguments
        class(lm_solver_options), intent(inout) :: x

        ! Set defaults
        x%method = FS_LEVENBERG_MARQUARDT_UPDATE
        x%finite_difference_step_size = sqrt(epsilon(1.0d0))
        x%damping_increase_factor = 11.0d0
        x%damping_decrease_factor = 9.0d0
    end subroutine

! ******************************************************************************
! OPERATORS
! ------------------------------------------------------------------------------
    module subroutine ic_equal(x, y)
        type(iteration_controls), intent(inout) :: x
        type(iteration_controls), intent(in) :: y
        x%max_iteration_count = y%max_iteration_count
        x%max_function_evaluations = y%max_function_evaluations
        x%gradient_tolerance = y%gradient_tolerance
        x%change_in_solution_tolerance = y%change_in_solution_tolerance
        x%residual_tolerance = y%residual_tolerance
        x%iteration_improvement_tolerance = y%iteration_improvement_tolerance
        x%max_iteration_between_updates = y%max_iteration_between_updates
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine ci_equal(x, y)
        type(convergence_info), intent(inout) :: x
        type(convergence_info), intent(in) :: y
        x%converge_on_gradient = y%converge_on_gradient
        x%gradient_value = y%gradient_value
        x%converge_on_solution_change = y%converge_on_solution_change
        x%solution_change_value = y%solution_change_value
        x%converge_on_residual_parameter = y%converge_on_residual_parameter
        x%residual_value = y%residual_value
        x%reach_iteration_limit = y%reach_iteration_limit
        x%iteration_count = y%iteration_count
        x%reach_function_evaluation_limit = y%reach_function_evaluation_limit
        x%function_evaluation_count = y%function_evaluation_count
        x%user_requested_stop = y%user_requested_stop
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine lso_equal(x, y)
        type(lm_solver_options), intent(inout) :: x
        type(lm_solver_options), intent(in) :: y
        x%method = y%method
        x%finite_difference_step_size = y%finite_difference_step_size
        x%damping_increase_factor = y%damping_increase_factor
        x%damping_decrease_factor = y%damping_decrease_factor
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
    ! Extracts the diagonal (D - MIN(M,N)) from a matrix (X - M-by-N).
    subroutine extract_diagonal(x, d)
        ! Arguments
        real(real64), intent(in) :: x(:,:)
        real(real64), intent(out) :: d(:)

        ! Local Variables
        integer(int32) :: i, m, n, mn

        ! Process
        m = size(x, 1)
        n = size(x, 2)
        mn = min(m, n)
        do i = 1, mn
            d(i) = x(i,i)
        end do
    end subroutine

! ------------------------------------------------------------------------------
end submodule
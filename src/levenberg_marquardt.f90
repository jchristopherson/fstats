submodule (fstats) levenberg_marquardt
! REFERENCES:
! 1. https://people.duke.edu/~hpgavin/ExperimentalSystems/lm.pdf
    use linalg
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
        integer(int32) :: m, n, flag, expected, actual, index
        real(real64), pointer :: f1p(:), f0p(:)
        real(real64), allocatable, target :: f1a(:), f0a(:), work(:)
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg
        
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
            go to 15
        end if
        if (present(f0)) then
            ! Check Size
            if (size(f0) /= m) then
                actual = size(f0)
                expected = m
                index = 2
                go to 10
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
                actual = size(f1)
                expected = m
                index = 3
                go to 10
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

        ! Array Size Error Handling
10      continue
        allocate(character(len = 512) :: errmsg)
        select case (index)
        case (1)
            write(errmsg, 100) "The dependent variable data array (", actual, &
                ") is not the same length as the independent data array (", &
                expected, ")."
        case (2)
            write(errmsg, 100) &
                "The function value array was expected to be of size ", &
                expected, ", but was found to be of size ", actual, "."
        case (3)
            write(errmsg, 100) &
                "The output function value array was expected to be of size ", &
                expected, ", but was found to be of size ", actual, "."
        end select
        call errmgr%report_error("regression_jacobian_1", trim(errmsg), &
            FS_ARRAY_SIZE_ERROR)
        return

        ! Jacobian Size Error Handling
15      continue
        allocate(character(len = 512) :: errmsg)
        write(errmsg, 101) "The Jacobian matrix was expected to be of size ", &
            m, "-by-", n, ", but was found to be ", size(jac, 1), "-by-", &
            size(jac, 2), "."
        call errmgr%report_error("regression_jacobian_1", trim(errmsg), &
            FS_ARRAY_SIZE_ERROR)
        return

        ! Memroy Allocation Error Handling
20      continue
        allocate(character(len = 512) :: errmsg)
        write(errmsg, 102) "Memory allocation error code ", flag, "."
        call errmgr%report_error("regression_jacobian_1", &
            trim(errmsg), ML_OUT_OF_MEMORY_ERROR)
        return

        ! Formatting
100     format(A, I0, A, I0, A)
101     format(A, I0, A, I0, A, I0, A, I0, A)
102     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------


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
        call recip_mult_array(h2, dy)   ! compute dy / h2
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
    subroutine lm_matrix(fun, xdata, ydata, pOld, yOld, dX2, jac, p, weights, &
        neval, update, step, JtWJ, JtWdy, X2, yNew, stop, work, mwork)
        ! Arguments
        procedure(regression_function), pointer :: fun
        real(real64), intent(in) :: xdata(:), ydata(:), pOld(:), yOld(:), &
            p(:), weights(:)
        real(real64), intent(in) :: dX2, step
        real(real64), intent(inout) :: jac(:,:)
        integer(int32), intent(inout) :: neval
        logical, intent(in) :: update
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
    ! - lambda:
    ! - maxP: maximum limits on the parameters.  Use huge() or larger for no constraints (N-by-1)
    ! - minP: minimum limits on the parameters.  Use -huge() or smaller for no constraints (N-by-1)
    ! - weights: A weighting vector (M-by-1)
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
        maxP, minP, weights, JtWJ, JtWdy, h, pNew, deltaY, yNew, X2, stop, &
        iwork, err)
        ! Arguments
        procedure(regression_function), pointer :: fun
        real(real64), intent(in) :: xdata(:), ydata(:), p(:), maxP(:), &
            minP(:), weights(:), JtWdy(:)
        real(real64), intent(in) :: lambda
        integer(int32), intent(inout) :: neval, niter
        integer(int32), intent(in) :: update
        real(real64), intent(inout) :: JtWJ(:,:)
        real(real64), intent(out) :: h(:), pNew(:), deltaY(:), yNew(:)
        real(real64), intent(out) :: X2
        logical, intent(out) :: stop
        integer(int32), intent(out) :: iwork(:)
        class(errors), intent(inout) :: err

        ! Local Variables
        integer(int32) :: i, n

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
    end subroutine

! ------------------------------------------------------------------------------
    ! A Levenberg-Marquardt solver.
    !
    ! Inputs:
    ! - fun: The function to evaluate
    ! - xdata: The independent coordinate data to fit (M-by-1)
    ! - ydata: The dependent coordinate data to fit (M-by-1)
    ! - p: current set of parameters (N-by-1)
    ! - weights:
    ! - maxP:
    ! - minP:
    ! - controls:
    ! - opt:
    !
    ! Outputs:
    ! - p:
    ! - y:
    ! - resid:
    ! - JtWJ:
    !
    ! - stop:
    ! - err: An error handling object
    subroutine lm_solve(fun, xdata, ydata, p, weights, maxP, minP, controls, opt, y, resid, JtWJ, info, stop, err)
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

        ! Local Variables
        logical :: update
        integer(int32) :: i, m, n, dof, flag, neval, niter
        real(real64) :: dX2, X2, X2Old, X2Try, lambda, alpha, dpJh, rho, nu, &
            step
        real(real64), allocatable :: pOld(:), yOld(:), J(:,:), JtWdy(:), &
            work(:), mwork(:,:), pTry(:), yTemp(:), JtWJc(:,:), h(:)
        integer(int32), allocatable :: iwork(:)
        character(len = :), allocatable :: errmsg

        ! Initialization
        update = .false.
        m = size(xdata)
        n = size(p)
        dof = m - n
        niter = 0
        step = opt%finite_difference_step_size
        stop = .false.

        ! Local Memory Allocation
        allocate(pOld(n), source = 0.0d0, stat = flag)
        if (flag == 0) allocate(yOld(m), source = 0.0d0, stat = flag)
        if (flag == 0) allocate(J(m, n), stat = flag)
        if (flag == 0) allocate(JtWdy(n), stat = flag)
        if (flag == 0) allocate(work(m + n), stat = flag)
        if (flag == 0) allocate(mwork(n, m), stat = flag)
        if (flag == 0) allocate(pTry(n), stat = flag)
        if (flag == 0) allocate(h(m), stat = flag)
        if (flag == 0) allocate(yTemp(m), stat = flag)
        if (flag == 0) allocate(JtWJc(n, n), stat = flag)
        if (flag == 0) allocate(iwork(n), stat = flag)
        if (flag /= 0) go to 10

        ! Perform an initial function evaluation
        call fun(xdata, p, y, stop)
        neval = 1

        ! Evaluate the problem matrices
        call lm_matrix(fun, xdata, ydata, pOld, yOld, 1.0d0, J, p, weights, &
            neval, .true., step, JtWJ, JtWdy, X2, y, stop, work, mwork)
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
                yTemp, X2Try, stop, iwork, err)
            if (stop) go to 5
            if (err%has_error_occurred()) return

            ! Perform a quadratic line update in the H direction
            if (opt%method == FS_QUADRATIC_UPDATE) then
                dpJh = dot_product(JtWdy, h)
                alpha = dpJh / (0.5d0 * (X2Try - X2) + 2.0d0 * dpJh)
                h = alpha * h

                do i = 1, n
                    pTry(i) = min(max(minP(i), p(i) + h(i)), maxP(i))
                end do

                call fun(xdata, pTry, yTemp, stop)
                if (stop) go to 5
                neval = neval + 1
                resid = ydata - yTemp
                X2Try = dot_product(resid, resid * weights)
            end if

            ! Compare the improvement in Chi-squared vs. the theoretical
            ! improvement of a single LM update
            if (opt%method == FS_LEVENBERG_MARQUARDT_UPDATE) then
                call extract_diagonal(JtWJ, work(1:n))
                work(1:n) = lambda * work(1:n) * h + JtWdy
            else
                work(1:n) = lambda * h + JtWdy
            end if
            rho = (X2 - X2Try) / abs(dot_product(h, work(1:n)))
            if (rho > controls%iteration_improvement_tolerance) then
                ! Things are getting better at an acceptable rate
                dX2 = X2 - X2Old
                X2Old = X2
                pOld = P
                yOld = y
                p = pTry

                ! Recompute the matrices
                call lm_matrix(fun, xdata, ydata, pOld, yOld, dX2, J, p, &
                    weights, neval, update, step, JtWJ, JtWdy, X2, y, stop, &
                    work, mwork)
                if (stop) go to 5
                JtWJc = JtWJ

                ! Decrease lambda
                select case (opt%method)
                case (FS_LEVENBERG_MARQUARDT_UPDATE)
                    lambda = max(lambda / opt%damping_decrease_factor, 1.0d-7)
                case (FS_QUADRATIC_UPDATE)
                    lambda = max(lambda / (1.0d0 + alpha), 1.0d-7)
                case (FS_NIELSEN_UPDATE)
                    lambda = lambda * max(1.0d0 / 3.0d0, &
                        1.0d0 - (2.0d0 * rho - 1.0d0)**3)
                    nu = 2.0d0
                end select
            else
                ! The iteration is not improving
                X2 = X2Old
                if (mod(niter, 2 * n) /= 0) then
                    ! Force a rank-1 update of the system matrices
                    call lm_matrix(fun, xdata, ydata, pOld, yOld, -1.0d0, J, &
                        p, weights, neval, update, step, JtWJ, JtWdy, dX2, &
                        y, stop, work, mwork)
                    if (stop) go to 5
                    JtWJc = JtWJ
                end if

                ! Increase lambda
                select case (opt%method)
                case (FS_LEVENBERG_MARQUARDT_UPDATE)
                    lambda = min(lambda * opt%damping_increase_factor, 1.0d7)
                case (FS_QUADRATIC_UPDATE)
                    lambda = lambda + abs((X2Try - X2) / 2.0d0 / alpha)
                case (FS_NIELSEN_UPDATE)
                    lambda = lambda * nu
                    nu = 2.0d0 * nu
                end select
            end if

            ! Determine the matrix update scheme
            update = mod(niter, 2 * n) > 0

            ! Test for convergence
            if (lm_check_convergence(controls, dof, resid, niter, neva, &
                JtWdy, h, p, X2, info)) &
            then
                exit main
            end if
        end do main

        ! End
        return

        ! User Requested End
5       continue
        return

        ! Memory Error Handling
10      continue
        allocate(character(len = 512) :: errmsg)
        write(errmsg, 100) "Memory allocation error code ", flag, "."
        call err%report_error("lm_solve", &
            trim(errmsg), ML_OUT_OF_MEMORY_ERROR)
        return

        ! Formatting
100     format(A, I0, A)
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
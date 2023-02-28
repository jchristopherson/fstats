submodule (fstats) levenberg_marquardt
! REFERENCES:
! 1. https://people.duke.edu/~hpgavin/ExperimentalSystems/lm.pdf
    use linalg
contains
! ------------------------------------------------------------------------------
    module subroutine regression_jacobian_1(fun, xdata, ydata, params, &
        jac, stop, f0, f1, step, err)
        ! Arguments
        procedure(regression_function), intent(in), pointer :: fun
        real(real64), intent(in) :: xdata(:), ydata(:), params(:)
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
        if (size(ydata) /= m) then
            actual = size(ydata)
            expected = m
            index = 1
            go to 10
        end if
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
            call fun(xdata, ydata, params, f0p, stop)
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
        call jacobian_finite_diff(fun, xdata, ydata, params, f0p, jac, f1p, &
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
    ! - ydata: The dependent coordinate data to fit (M-by-1)
    ! - params: The model parameters (N-by-1)
    ! - f0: The current model estimate (M-by-1)
    ! - step: The differentiation step size
    !
    ! Outputs:
    ! - jac: The Jacobian matrix (M-by-N)
    ! - f1: A workspace array for the model output (M-by-1)
    ! - stop: A flag allowing the user to terminate model execution
    ! - work: A workspace array for the model parameters (N-by-1)
    subroutine jacobian_finite_diff(fun, xdata, ydata, params, f0, jac, f1, &
        stop, step, work)
        ! Arguments
        procedure(regression_function), intent(in), pointer :: fun
        real(real64), intent(in) :: xdata(:), ydata(:), params(:)
        real(real64), intent(inout) :: f0(:)
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
            call fun(xdata, ydata, work, f1, stop)
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
    ! - weights: A weighting vector (N-by-1)
    ! - neval: Current number of function evaluations
    !
    ! Outputs:
    ! - JtWJ: linearized Hessian matrix (inverse of the covariance matrix) (N-by-N)
    ! - JtWdy: linearized fitting vector (N-by-M)
    ! - X2: Updated Chi-squared criteria
    ! - yNew: model evaluated with parameters of p (M-by-1)
    ! - jac: updated Jacobian matrix (M-by-N)
    ! - neval: updated count of function evaluations
    subroutine lm_matrix(fun, xdata, ydata, pOld, yOld, dX2, jac, p, weights, neval, JtwJ, JtWdy, X2, yNew)
        ! Arguments
        procedure(regression_function), pointer :: fun
        real(real64), intent(in) :: xdata(:), ydata(:), pOld(:), yOld(:), p(:), weights(:)
        real(real64), intent(in) :: dX2
        real(real64), intent(inout) :: jac(:,:)
        integer(int32), intent(inout) :: neval
        real(real64), intent(out) :: JtWJ(:,:), JtWdy(:,:)
        real(real64), intent(out) :: X2

        ! Local Variables

        ! Process
    end subroutine

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
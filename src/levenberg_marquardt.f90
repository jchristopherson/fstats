submodule (fstats) levenberg_marquardt
! REFERENCES:
! 1. https://people.duke.edu/~hpgavin/ExperimentalSystems/lm.pdf
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


! ******************************************************************************
! PRIVATE ROUTINES
! ------------------------------------------------------------------------------
    subroutine jacobian_finite_diff(fun, xdata, ydata, params, f0, jac, f1, stop, step, work)
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
end submodule
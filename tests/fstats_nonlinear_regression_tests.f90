module fstats_nonlinear_regression_tests
    use iso_fortran_env
    use fstats_test_helper
    use fstats
    implicit none

    interface
        function prototype(x) result(y)
            use iso_fortran_env, only : real64
            real(real64), intent(in) :: x
            real(real64) :: y
        end function
    end interface
contains
! ------------------------------------------------------------------------------
    function square_x(x) result(y)
        real(real64), intent(in) :: x
        real(real64) :: y
        y = x**2
    end function

    function test_prototype_function_call() result(rst)
        ! Arguments
        logical :: rst

        ! Local Variables
        real(real64), parameter :: x = 4.0d0
        real(real64), parameter :: ans = x**2
        procedure(prototype), pointer :: fun
        real(real64) :: y

        ! Initialization
        rst = .true.
        fun => square_x

        ! Process
        y = fun(x)
        if (.not.is_equal(y, ans)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_prototype_function_call 1-1"
        end if
    end function

! ------------------------------------------------------------------------------
    subroutine jac_fcn_1(xdata, params, f, stop)
        ! Arguments
        real(real64), intent(in) :: xdata(:), params(:)
        real(real64), intent(out) :: f(:)
        logical, intent(out) :: stop

        ! Function
        ! f(x) = a1 * x**2 + a2 * x + a3
        f = params(1) * xdata**2 + params(2) * xdata + params(3)

        ! End
        stop = .false.
    end subroutine

    function test_jacobian() result(rst)
        ! Arguments
        logical :: rst

        ! Local Variables
        real(real64), parameter :: tol = 1.0d-4
        real(real64), parameter :: a1 = 0.5d0
        real(real64), parameter :: a2 = -1.35d0
        real(real64), parameter :: a3 = 1.2d0
        real(real64) :: x(3), ans(3,3), jac(3,3), params(3)
        procedure(regression_function), pointer :: fun
        logical :: stop

        ! Initialization
        rst = .true.
        x = [-1.25d0, 0.0d0, 1.25d0]
        params = [a1, a2, a3]

        ! Evaluate the real Jacobian: J(i,j) = dF(i) / dx(j)
        ans = reshape([ &
            x(1)**2, x(2)**2, x(3)**2, &
            x(1), x(2), x(3), &
            1.0d0, 1.0d0, 1.0d0], &
            [3, 3] &
        )

        ! Compute the Jacobian numerically
        fun => jac_fcn_1
        call jacobian(fun, x, params, jac, stop)

        ! Test
        if (.not.is_equal(ans, jac, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_jacobian 1-1"
        end if
    end function

! ------------------------------------------------------------------------------
end module
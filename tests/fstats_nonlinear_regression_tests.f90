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
    subroutine nl_test_fun(x, p, f, stop)
        ! Arguments
        real(real64), intent(in) :: x(:), p(:)
        real(real64), intent(out) :: f(:)
        logical, intent(out) :: stop

        ! Function
        f = p(1) * x**3 + p(2) * x**2 + p(3) * x + p(4)

        ! End
        stop = .false.
    end subroutine

    function test_nl_least_squares() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: tol = 1.0d-5

        ! Local Variables
        procedure(regression_function), pointer :: fun
        real(real64) :: xp(21), yp(21), params(4), ymod(21), resid(21), ans(4)
        type(lm_solver_options) :: opt
        type(iteration_controls) :: tols

        ! Data to fit
        xp = [0.0d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.8d0, &
            0.9d0, 1.0d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0, 1.5d0, 1.6d0, 1.7d0, &
            1.8d0, 1.9d0, 2.0d0]
        yp = [1.216737514d0, 1.250032542d0, 1.305579195d0, 1.040182335d0, &
            1.751867738d0, 1.109716707d0, 2.018141531d0, 1.992418729d0, &
            1.807916923d0, 2.078806005d0, 2.698801324d0, 2.644662712d0, &
            3.412756702d0, 4.406137221d0, 4.567156645d0, 4.999550779d0, &
            5.652854194d0, 6.784320119d0, 8.307936836d0, 8.395126494d0, &
            10.30252404d0]

        ! Expected Solution:
        ans = [1.0647627571d0, -0.1223202909d0, 0.4466134462d0, 1.1866142244d0]

        ! Initialization
        fun => nl_test_fun
        params = 1.0d0

        ! Solve the problem
        call nonlinear_least_squares(fun, xp, yp, params, ymod, resid)

        ! Check the result
        if (is_equal(params, ans, tol)) then
            rst = .true.
        else
            rst = .false.
            print '(A)', "TEST FAILED: test_nl_least_squares 1-1"
        end if

        ! Test with a different options set
        call opt%set_to_default()
        opt%method = FS_NIELSEN_UPDATE
        params = 1.0d0
        call nonlinear_least_squares(fun, xp, yp, params, ymod, resid, &
            settings = opt)

        ! Check the result
        if (.not.is_equal(params, ans, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_nl_least_squares 1-2"
        end if

        ! Now test the quadratic update option
        call tols%set_to_default()
        tols%change_in_solution_tolerance = 1.0e-11
        opt%method = FS_QUADRATIC_UPDATE
        params = 1.0d0
        call nonlinear_least_squares(fun, xp, yp, params, ymod, resid, &
            settings = opt, controls = tols)

        ! Check the result
        if (.not.is_equal(params, ans, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_nl_least_squares 1-3"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_bootstrap_nl_least_squares() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: tol = 1.0d-6

        ! Local Variables
        procedure(regression_function), pointer :: fun
        real(real64) :: xp(21), yp(21), params(4), ymod(21), resid(21), ans(4)
        type(bootstrap_regression_statistics) :: stats(4)

        ! Data to fit
        xp = [0.0d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.8d0, &
            0.9d0, 1.0d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0, 1.5d0, 1.6d0, 1.7d0, &
            1.8d0, 1.9d0, 2.0d0]
        yp = [1.216737514d0, 1.250032542d0, 1.305579195d0, 1.040182335d0, &
            1.751867738d0, 1.109716707d0, 2.018141531d0, 1.992418729d0, &
            1.807916923d0, 2.078806005d0, 2.698801324d0, 2.644662712d0, &
            3.412756702d0, 4.406137221d0, 4.567156645d0, 4.999550779d0, &
            5.652854194d0, 6.784320119d0, 8.307936836d0, 8.395126494d0, &
            10.30252404d0]

        ! Expected Solution:
        ans = [1.0647627571d0, -0.1223202909d0, 0.4466134462d0, 1.1866142244d0]

        ! Initialization
        rst = .true.
        fun => nl_test_fun
        params = 1.0d0

        ! Solve the problem
        call bootstrap_nonlinear_least_squares(fun, xp, yp, params, ymod, &
            resid, stats = stats)
        
        ! Check the result
        if (.not.is_equal(params, ans, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_bootstrap_nl_least_squares 1-1"
        end if
    end function

! ------------------------------------------------------------------------------
end module
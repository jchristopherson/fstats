module fstats_interp_tests
    use iso_fortran_env
    use fortran_test_helper
    use fstats
    implicit none

contains
! ------------------------------------------------------------------------------
pure function linear_interp_reference(x, x1, x2, y1, y2) result(rst)
    real(real64), intent(in) :: x, x1, x2, y1, y2
    real(real64) :: rst

    rst = -((x1 - x) * y2 + (x - x2) * y1) / (x2 - x1)
end function

! ------------------------------------------------------------------------------
function test_linear_interp() result(rst)
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-8
    real(real64), parameter :: x(4) = [-1.5d0, 0.0d0, 1.5d0, 3.0d0]
    real(real64), parameter :: x1 = -0.2d0
    real(real64), parameter :: x2 = 2.75d0
    real(real64), parameter :: x3 = -2.0d0
    real(real64), parameter :: x4 = 4.0d0

    ! Local Variables
    real(real64) :: y(4), y1, y2, y3, y4, a1, a2, a3, a4
    type(linear_interpolator) :: interp

    ! Initialization
    rst = .true.
    call random_number(y)
    call interp%initialize(x, y)

    ! Test 1
    a1 = linear_interp_reference(x1, x(1), x(2), y(1), y(2))
    y1 = interp%interpolate_value(x1)
    if (.not.assert(a1, y1, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_linear_interp -1"
    end if

    ! Test 2
    a2 = linear_interp_reference(x2, x(3), x(4), y(3), y(4))
    y2 = interp%interpolate_value(x2)
    if (.not.assert(a2, y2, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_linear_interp -2"
    end if

    ! Test 3
    a3 = linear_interp_reference(x3, x(1), x(2), y(1), y(2))
    y3 = interp%interpolate_value(x3)
    if (.not.assert(a3, y3, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_linear_interp -3"
    end if

    ! Test 4
    a4 = linear_interp_reference(x4, x(3), x(4), y(3), y(4))
    y4 = interp%interpolate_value(x4)
    if (.not.assert(a4, y4, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_linear_interp -4"
    end if
end function

! ------------------------------------------------------------------------------
pure function quadratic_interp_reference(x, x1, x2, x3, y1, y2, y3) result(rst)
    real(real64), intent(in) :: x, x1, x2, x3, y1, y2, y3
    real(real64) :: rst

    real(real64) :: denom, a, b, c

    denom = (x2 - x1) * (x3 - x1) * (x3 - x2)
    a = (x2 * (y3 - y1) + x1 * (y2 - y3) + x3 * (y1 - y2)) / denom
    b = -(x2**2 * (y3 - y1) + x1**2 * (y2 - y3) + x3**2 * (y1 - y2)) / denom
    c = (y3 * (x1 * x2**2 - x1**2 * x2) + y2 * (x1**2 * x3 - x1 * x3**2) + &
        y1 * (x2 * x3**2 - x2**2 * x3)) / denom
    rst = a * x**2 + b * x + c
end function

! ------------------------------------------------------------------------------
function test_polynomial_interp() result(rst)
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-8
    real(real64), parameter :: x(4) = [-1.5d0, 0.0d0, 1.5d0, 3.0d0]
    real(real64), parameter :: x1 = -0.2d0
    real(real64), parameter :: x2 = 2.75d0
    real(real64), parameter :: x3 = -2.0d0
    real(real64), parameter :: x4 = 4.0d0

    ! Local Variables
    real(real64) :: y(4), y1, y2, y3, y4, a1, a2, a3, a4
    type(polynomial_interpolator) :: interp

    ! Initialization
    rst = .true.
    call random_number(y)
    call interp%initialize(2, x, y)

    ! Test 1
    a1 = quadratic_interp_reference(x1, x(1), x(2), x(3), y(1), y(2), y(3))
    y1 = interp%interpolate_value(x1)
    if (.not.assert(a1, y1, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_polynomial_interp -1"
    end if

    ! Test 2
    a2 = quadratic_interp_reference(x2, x(2), x(3), x(4), y(2), y(3), y(4))
    y2 = interp%interpolate_value(x2)
    if (.not.assert(a2, y2, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_polynomial_interp -2"
    end if

    ! Test 3
    a3 = quadratic_interp_reference(x3, x(1), x(2), x(3), y(1), y(2), y(3))
    y3 = interp%interpolate_value(x3)
    if (.not.assert(a3, y3, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_polynomial_interp -3"
    end if

    ! Test 4
    a4 = quadratic_interp_reference(x4, x(2), x(3), x(4), y(2), y(3), y(4))
    y4 = interp%interpolate_value(x4)
    if (.not.assert(a4, y4, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_polynomial_interp -4"
    end if
end function

! ------------------------------------------------------------------------------
pure function ref_poly(x, k) result(rst)
    !! Evaluates 1 + 2 x + 3 x**2 + ... + (k + 1) x**k.
    real(real64), intent(in) :: x
    integer(int32), intent(in) :: k
    real(real64) :: rst

    integer(int32) :: i

    rst = 0.0d0
    do i = 0, k
        rst = rst + real(i + 1, real64) * x**i
    end do
end function

! ------------------------------------------------------------------------------
pure function ref_poly_derivative(x, k) result(rst)
    !! Evaluates the derivative of ref_poly.
    real(real64), intent(in) :: x
    integer(int32), intent(in) :: k
    real(real64) :: rst

    integer(int32) :: i

    rst = 0.0d0
    do i = 1, k
        rst = rst + real(i, real64) * real(i + 1, real64) * x**(i-1)
    end do
end function

! ------------------------------------------------------------------------------
function test_polynomial_interp_orders() result(rst)
    !! An order K interpolant must reproduce a polynomial of degree K exactly,
    !! and must do so over the whole data range rather than just the interior.
    logical :: rst

    ! Parameters
    integer(int32), parameter :: n = 12
    real(real64), parameter :: tol = 1.0d-9

    ! Local Variables
    integer(int32) :: i, k
    real(real64) :: x(n), y(n), xi, ans
    type(polynomial_interpolator) :: interp

    ! Initialization
    rst = .true.

    do k = 1, 5
        do i = 1, n
            x(i) = real(i - 1, real64)
            y(i) = ref_poly(x(i), k)
        end do
        call interp%initialize(k, x, y)
        do i = 0, 220
            xi = 5.0d-2 * i
            ans = ref_poly(xi, k)
            if (abs(interp%interpolate_value(xi) - ans) > tol * abs(ans)) then
                rst = .false.
                print "(A, I0, A, F8.3)", &
                    "TEST FAILED: test_polynomial_interp_orders, order ", k, &
                    ", x = ", xi
                exit
            end if
        end do
    end do

    ! An order 1 interpolant must use the interval that brackets the point,
    ! not a neighboring one
    do i = 1, n
        x(i) = real(i - 1, real64)
        y(i) = x(i)**2
    end do
    call interp%initialize(1, x, y)
    do i = 1, n - 1
        xi = x(i) + 5.0d-1
        ans = y(i) + (xi - x(i)) * (y(i+1) - y(i)) / (x(i+1) - x(i))
        if (.not.assert(interp%interpolate_value(xi), ans, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_polynomial_interp_orders -2"
            exit
        end if
    end do
end function

! ------------------------------------------------------------------------------
function test_spline_interp() result(rst)
    logical :: rst

    ! Parameters
    integer(int32), parameter :: n = 11
    real(real64), parameter :: tol = 1.0d-10

    ! Local Variables
    integer(int32) :: i
    real(real64) :: x(n), y(n), xi, ans
    type(spline_interpolator) :: interp

    ! Initialization
    rst = .true.
    do i = 1, n
        x(i) = real(i - 1, real64)
        y(i) = ref_poly(x(i), 3)
    end do

    ! Test 1 - a cubic spline with the exact end slopes reproduces a cubic
    call interp%initialize(x, y, SPLINE_KNOWN_FIRST_DERIVATIVE, &
        ref_poly_derivative(x(1), 3), SPLINE_KNOWN_FIRST_DERIVATIVE, &
        ref_poly_derivative(x(n), 3))
    do i = 0, 200
        xi = 5.0d-2 * i
        ans = ref_poly(xi, 3)
        if (abs(interp%interpolate_value(xi) - ans) > tol * abs(ans)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_spline_interp -1"
            exit
        end if
    end do

    ! Test 2 - the not-a-knot condition also reproduces a cubic, and exercises
    ! the genuinely pentadiagonal path
    call interp%initialize(x, y, SPLINE_CONTINUOUS_THIRD_DERIVATIVE, 0.0d0, &
        SPLINE_CONTINUOUS_THIRD_DERIVATIVE, 0.0d0)
    do i = 0, 200
        xi = 5.0d-2 * i
        ans = ref_poly(xi, 3)
        if (abs(interp%interpolate_value(xi) - ans) > tol * abs(ans)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_spline_interp -2"
            exit
        end if
    end do

    ! Test 3 - the default end condition is quadratic over the end intervals,
    ! so a quadratic must be reproduced
    do i = 1, n
        y(i) = ref_poly(x(i), 2)
    end do
    call interp%initialize(x, y)
    do i = 0, 200
        xi = 5.0d-2 * i
        ans = ref_poly(xi, 2)
        if (abs(interp%interpolate_value(xi) - ans) > tol * abs(ans)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_spline_interp -3"
            exit
        end if
    end do

    ! Test 4 - the nodes must be reproduced
    do i = 1, n
        if (.not.assert(interp%interpolate_value(x(i)), y(i), tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_spline_interp -4"
            exit
        end if
    end do
end function

! ------------------------------------------------------------------------------
function test_hermite_interp() result(rst)
    !! Given N points the Hermite interpolant is of degree 2N - 1, so with
    !! N = 4 a degree 7 polynomial and its derivative must be reproduced.
    logical :: rst

    ! Parameters
    integer(int32), parameter :: n = 4
    integer(int32), parameter :: nq = 31
    real(real64), parameter :: tol = 1.0d-10

    ! Local Variables
    integer(int32) :: i
    real(real64) :: x(n), y(n), yp(n), xq(nq), v(nq), d(nq), ans
    type(hermite_interpolator) :: interp

    ! Initialization
    rst = .true.
    do i = 1, n
        x(i) = real(i - 1, real64)
        y(i) = ref_poly(x(i), 7)
        yp(i) = ref_poly_derivative(x(i), 7)
    end do
    call interp%initialize(x, y, yp)

    do i = 1, nq
        xq(i) = 1.0d-1 * (i - 1)
    end do
    call interp%interpolate_with_derivative(xq, v, d)

    ! Test 1 - the values
    do i = 1, nq
        ans = ref_poly(xq(i), 7)
        if (abs(v(i) - ans) > tol * abs(ans)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_hermite_interp -1"
            exit
        end if
    end do

    ! Test 2 - the derivatives
    do i = 1, nq
        ans = ref_poly_derivative(xq(i), 7)
        if (abs(d(i) - ans) > tol * abs(ans)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_hermite_interp -2"
            exit
        end if
    end do

    ! Test 3 - the node values and derivatives must be matched exactly
    block
        real(real64) :: vn(n), dn(n)
        call interp%interpolate_with_derivative(x, vn, dn)
        do i = 1, n
            if (.not.assert(vn(i), y(i), tol) .or. &
                .not.assert(dn(i), yp(i), tol)) then
                rst = .false.
                print "(A)", "TEST FAILED: test_hermite_interp -3"
                exit
            end if
        end do
    end block
end function

! ------------------------------------------------------------------------------
end module
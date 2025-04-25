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

! ------------------------------------------------------------------------------
end module
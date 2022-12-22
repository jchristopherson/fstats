module fstats_experimental_design_tests
    use iso_fortran_env
    use fstats
    use fstats_test_helper
    implicit none
contains
! ------------------------------------------------------------------------------
    function get_full_matrix_size_test_1() result(rst)
        ! Arguments
        logical :: rst

        ! Local Variables
        integer(int32), parameter :: rowcount = 24
        integer(int32), parameter :: columncount = 3
        integer(int32) :: vars(columncount), m, n

        ! Process
        rst = .true.
        vars = [2, 4, 3]
        call get_full_factorial_matrix_size(vars, m, n)

        if (m /= rowcount) then
            rst = .false.
            print '(A)', "TEST FAILED: get_full_matrix_size_test 1 - 1"
        end if

        if (n /= columncount) then
            rst = .false.
            print '(A)', "TEST FAILED: get_full_matrix_size_test 1 - 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function full_factorial_test_1() result(rst)
        ! Arguments
        logical :: rst

        ! Local Variables
        integer(int32) :: vars(3), tbl(24, 3), ans(24, 3)

        ! Initialization
        rst = .true.
        vars = [2, 4, 3]
        ans = reshape([ &
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
            1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, &
            1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3], &
            [24, 3] &
        )

        ! Process
        call full_factorial(vars, tbl)
        if (.not.is_equal(tbl, ans)) then
            rst = .false.
            print '(A)', "TEST FAILED: Full Factorial Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
end module
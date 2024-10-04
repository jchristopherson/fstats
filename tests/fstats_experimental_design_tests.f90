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
    function test_eval_model_main_effects() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        integer(int32), parameter :: nfactor = 3
        integer(int32), parameter :: npts = 50
        integer(int32), parameter :: nparams = 1 + nfactor

        ! Local Variables
        integer(int32) :: i
        real(real64) :: x(npts, nfactor), mdl(nparams), y(npts), ans(npts)

        ! Initialization
        rst = .true.
        call random_number(mdl)
        call random_number(x)
        
        ! Define the answer
        ans = mdl(1)
        do i = 1, nfactor
            ans = ans + mdl(i+1) * x(:,i)
        end do

        ! Perform the calculation
        y = doe_evaluate_model(1, mdl, x)

        ! Test
        if (.not.is_equal(y, ans)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_eval_model_main_effects -1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_eval_model_two_way() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        integer(int32), parameter :: nfactor = 3
        integer(int32), parameter :: npts = 50
        integer(int32), parameter :: nparams = 1 + nfactor + nfactor * (nfactor - 1)

        ! Local Variables
        real(real64) :: x(npts,nfactor), mdl(nparams), y(npts), ans(npts)

        ! Initialization
        rst = .true.
        call random_number(x)
        call random_number(mdl)

        ! Define the answer
        ans = mdl(1) + mdl(2) * x(:,1) + mdl(3) * x(:,2) + mdl(4) * x(:,3) + &
            mdl(5) * x(:,1) * x(:,2) + mdl(6) * x(:,1) * x(:,3) + &
            mdl(7) * x(:,2) * x(:,1) + mdl(8) * x(:,2) * x(:,3) + &
            mdl(9) * x(:,3) * x(:,1) + mdl(10) * x(:,3) * x(:,2)

        ! Perform the calculation
        y = doe_evaluate_model(2, mdl, x)

        ! Test
        if (.not.is_equal(y, ans)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_eval_model_two_way -1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_eval_model_three_way() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        integer(int32), parameter :: nfactor = 3
        integer(int32), parameter :: npts = 50
        integer(int32), parameter :: nparams = 1 + nfactor + &
            nfactor * (nfactor - 1) + nfactor * (nfactor**2 - 1)

        ! Local Variables
        real(real64) :: x(npts,nfactor), mdl(nparams), y(npts), ans(npts)

        ! Initialization
        rst = .true.
        call random_number(x)
        call random_number(mdl)

        ! Define the answer
        ans = mdl(1) + mdl(2) * x(:,1) + mdl(3) * x(:,2) + mdl(4) * x(:,3) + &
            mdl(5) * x(:,1) * x(:,2) + mdl(6) * x(:,1) * x(:,3) + &
            mdl(7) * x(:,2) * x(:,1) + mdl(8) * x(:,2) * x(:,3) + &
            mdl(9) * x(:,3) * x(:,1) + mdl(10) * x(:,3) * x(:,2) + &
            mdl(11) * x(:,1) * x(:,1) * x(:,2) + &
            mdl(12) * x(:,1) * x(:,1) * x(:,3) + &
            mdl(13) * x(:,1) * x(:,2) * x(:,1) + &
            mdl(14) * x(:,1) * x(:,2) * x(:,2) + &
            mdl(15) * x(:,1) * x(:,2) * x(:,3) + &
            mdl(16) * x(:,1) * x(:,3) * x(:,1) + &
            mdl(17) * x(:,1) * x(:,3) * x(:,2) + &
            mdl(18) * x(:,1) * x(:,3) * x(:,3) + &
            mdl(19) * x(:,2) * x(:,1) * x(:,1) + &
            mdl(20) * x(:,2) * x(:,1) * x(:,2) + &
            mdl(21) * x(:,2) * x(:,1) * x(:,3) + &
            mdl(22) * x(:,2) * x(:,2) * x(:,1) + &
            mdl(23) * x(:,2) * x(:,2) * x(:,3) + &
            mdl(24) * x(:,2) * x(:,3) * x(:,1) + &
            mdl(25) * x(:,2) * x(:,3) * x(:,2) + &
            mdl(26) * x(:,2) * x(:,3) * x(:,3) + &
            mdl(27) * x(:,3) * x(:,1) * x(:,1) + &
            mdl(28) * x(:,3) * x(:,1) * x(:,2) + &
            mdl(29) * x(:,3) * x(:,1) * x(:,3) + &
            mdl(30) * x(:,3) * x(:,2) * x(:,1) + &
            mdl(31) * x(:,3) * x(:,2) * x(:,2) + &
            mdl(32) * x(:,3) * x(:,2) * x(:,3) + &
            mdl(33) * x(:,3) * x(:,3) * x(:,1) + &
            mdl(34) * x(:,3) * x(:,3) * x(:,2)


        ! Perform the calculation
        y = doe_evaluate_model(3, mdl, x)

        ! Test
        if (.not.is_equal(y, ans)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_eval_model_three_way -1"
        end if
    end function

! ------------------------------------------------------------------------------
end module
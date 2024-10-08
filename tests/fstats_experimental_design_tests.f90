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
    function test_doe_fit_1() result(rst)
        ! Arguments
        logical :: rst

        ! Variables
        real(real64), parameter :: tol = 1.0d-6
        integer(int32), parameter :: npts = 31
        integer(int32) :: i
        real(real64) :: x(npts), y(npts), ymod(npts), resid(npts), lsmdl(2), &
            xd(npts,1)
        type(regression_statistics) :: lsstats(2)
        type(doe_model) :: mdl

        ! Initialization
        rst = .true.
        x = [ &
            0.0d0, &
            0.1d0, &
            0.2d0, &
            0.3d0, &
            0.4d0, &
            0.5d0, &
            0.6d0, &
            0.7d0, &
            0.8d0, &
            0.9d0, &
            1.0d0, &
            1.1d0, &
            1.2d0, &
            1.3d0, &
            1.4d0, &
            1.5d0, &
            1.6d0, &
            1.7d0, &
            1.8d0, &
            1.9d0, &
            2.0d0, &
            2.1d0, &
            2.2d0, &
            2.3d0, &
            2.4d0, &
            2.5d0, &
            2.6d0, &
            2.7d0, &
            2.8d0, &
            2.9d0, &
            3.0d0 &
        ]
        y = [ &
            0.577855138260449d0, &
            0.614883095604222d0, &
            0.633891127488559d0, &
            0.718405829701721d0, &
            0.753668502759107d0, &
            0.814967857310145d0, &
            0.861870996499704d0, &
            0.925100533744381d0, &
            0.947038018520063d0, &
            1.025198043343280d0, &
            1.042142354497610d0, &
            1.121528566784440d0, &
            1.177570314994070d0, &
            1.229237567525370d0, &
            1.261114062593870d0, &
            1.296408162551430d0, &
            1.394353657051120d0, &
            1.367144391560370d0, &
            1.428164431435150d0, &
            1.548944935073270d0, &
            1.505100149282990d0, &
            1.560701023751520d0, &
            1.609113012481530d0, &
            1.663687366875500d0, &
            1.707149545456870d0, &
            1.800935947618110d0, &
            1.847819988906440d0, &
            1.884242821675810d0, &
            1.966174239373140d0, &
            1.977005266443110d0, &
            2.034137257154140d0 &    
        ]
        xd(:,1) = x

        ! Use a straight least-squares method to fit the data
        call linear_least_squares(1, .true., x, y, lsmdl, ymod, resid, &
            stats = lsstats)

        ! Fit the model using the DOE routine
        mdl = doe_fit_model(1, xd, y)

        ! Compare coefficients
        if (.not.is_equal(lsmdl, mdl%coefficients, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_doe_fit_1 -1"
        end if

        ! Compare statistics
        do i = 1, size(lsstats)
            if (.not.is_equal(lsstats(i)%confidence_interval, &
                mdl%stats(i)%confidence_interval, tol)) &
            then
                rst = .false.
                print 100, "TEST FAILED: test_doe_fit_1 -2 ", i
            end if

            if (.not.is_equal(lsstats(i)%probability, &
                mdl%stats(i)%probability, tol)) &
            then
                rst = .false.
                print 100, "TEST FAILED: test_doe_fit_1 -3 ", i
            end if

            if (.not.is_equal(lsstats(i)%standard_error, &
                mdl%stats(i)%standard_error, tol)) &
            then
                rst = .false.
                print 100, "TEST FAILED: test_doe_fit_1 -4 ", i
            end if

            if (.not.is_equal(lsstats(i)%t_statistic, &
                mdl%stats(i)%t_statistic, tol)) &
            then
                rst = .false.
                print 100, "TEST FAILED: test_doe_fit_1 -5 ", i
            end if
        end do

        ! Formatting
    100 format(A, I0)
    end function

! ------------------------------------------------------------------------------
end module
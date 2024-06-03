module fstats_regression_tests
    use iso_fortran_env
    use fstats
    use fstats_test_helper
    implicit none
contains
! ------------------------------------------------------------------------------
    function design_matrix_test_1() result(rst)
        ! Arguments
        logical :: rst

        ! Variables
        integer(int32), parameter :: npts = 50
        integer(int32), parameter :: order1 = 1
        integer(int32), parameter :: order2 = 4
        real(real64), parameter :: dx = 1.0d-2
        integer(int32) :: i
        real(real64) :: x(npts), c1(npts, 2), c2(npts, 1), c3(npts, 5), &
            ans1(npts, 2), ans2(npts, 1), ans3(npts, 5)

        ! Initialization
        rst = .true.
        x(1) = 0.0d0
        do i = 2, npts
            x(i) = x(i-1) + dx
        end do

        ans1(:,1) = 1.0d0
        ans1(:,2) = x

        ans2(:,1) = x

        ans3(:,1) = 1.0d0
        ans3(:,2) = x
        ans3(:,3) = x**2
        ans3(:,4) = x**3
        ans3(:,5) = x**4

        ! Test 1 - linear w/ intercept
        call design_matrix(order1, .true., x, c1)
        if (.not.is_equal(c1, ans1)) then
            rst = .false.
            print '(A)', "TEST FAILED: Design Matrix Test 1 - 1"
        end if

        ! Test 2 - linear w/o intercept
        call design_matrix(order1, .false., x, c2)
        if (.not.is_equal(c2, ans2)) then
            rst = .false.
            print '(A)', "TEST FAILED: Design Matrix Test 1 - 2"
        end if

        ! Test 3 - 4th order w/ intercept
        call design_matrix(order2, .true., x, c3)
        if (.not.is_equal(c3, ans3)) then
            rst = .false.
            print '(A)', "TEST FAILED: Design Matrix Test 1 - 3"
        end if
    end function

! ------------------------------------------------------------------------------
    function regression_test_1() result(rst)
        ! Arguments
        logical :: rst

        ! Variables
        integer(int32), parameter :: npts = 31
        real(real64), parameter :: tol = 1.0d-3
        real(real64), parameter :: se1 = 8.0d-3
        real(real64), parameter :: se2 = 5.0d-3
        real(real64), parameter :: t1 = 67.836d0
        real(real64), parameter :: t2 = 101.127d0
        real(real64), parameter :: slope2 = 0.7668972622d0
        real(real64) :: x(npts), y(npts), ymod(npts), resid(npts), c1(2), c2(2)
        type(regression_statistics) :: s1(2), s2(1)

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

        ! Fit the model - linear model with intercept
        call linear_least_squares(1, .true., x, y, c1, ymod, resid, s1)

        if (.not.is_equal(s1(1)%standard_error, se1, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: Regression Test 1 - 1"
        end if
        if (.not.is_equal(s1(2)%standard_error, se2, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: Regression Test 1 - 2"
        end if
        if (.not.is_equal(s1(1)%t_statistic, t1, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: Regression Test 1 - 3"
        end if
        if (.not.is_equal(s1(2)%t_statistic, t2, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: Regression Test 1 - 4"
        end if

        ! Fit the model - linear model without intercept
        call linear_least_squares(1, .false., x, y, c2, ymod, resid, s2)

        if (.not.is_equal(c2(1), 0.0d0)) then
            rst = .false.
            print '(A)', "TEST FAILED: Regression Test 1 - 5"
        end if
        if (.not.is_equal(c2(2), slope2)) then
            rst = .false.
            print '(A)', "TEST FAILED: Regression Test 1 - 6"
        end if
    end function

! ------------------------------------------------------------------------------
end module
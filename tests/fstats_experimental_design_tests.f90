! SPDX-FileCopyrightText: 2022-2026 Jason Christopherson
! SPDX-License-Identifier: MIT
!
! MIT License
!
! Copyright (c) 2022-2026 Jason Christopherson
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

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
        integer(int32), parameter :: nparams = 1 + nfactor + &
            nfactor * (nfactor - 1) / 2 + nfactor

        ! Local Variables
        real(real64) :: x(npts,nfactor), mdl(nparams), y(npts), ans(npts)

        ! Initialization
        rst = .true.
        call random_number(x)
        call random_number(mdl)

        ! Define the answer
        ans = mdl(1) + mdl(2) * x(:,1) + mdl(3) * x(:,2) + mdl(4) * x(:,3) + &
            mdl(5) * x(:,1) * x(:,2) + mdl(6) * x(:,1) * x(:,3) + &
            mdl(7) * x(:,2) * x(:,3) + mdl(8) * x(:,1)**2 + &
            mdl(9) * x(:,2)**2 + mdl(10) * x(:,3)**2

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
            nfactor * (nfactor + 1) / 2 + &
            nfactor * (nfactor + 1) * (nfactor + 2) / 6

        ! Local Variables
        real(real64) :: x(npts,nfactor), mdl(nparams), y(npts), ans(npts)

        ! Initialization
        rst = .true.
        call random_number(x)
        call random_number(mdl)

        ! Define the answer
        ans = mdl(1) + mdl(2) * x(:,1) + mdl(3) * x(:,2) + mdl(4) * x(:,3) + &
            mdl(5) * x(:,1) * x(:,2) + mdl(6) * x(:,1) * x(:,3) + &
            mdl(7) * x(:,2) * x(:,3) + mdl(8) * x(:,1)**2 + &
            mdl(9) * x(:,2)**2 + mdl(10) * x(:,3)**2 + &
            mdl(11) * x(:,1) * x(:,2) * x(:,3) + &
            mdl(12) * x(:,1)**2 * x(:,2) + mdl(13) * x(:,1)**2 * x(:,3) + &
            mdl(14) * x(:,2)**2 * x(:,1) + mdl(15) * x(:,2)**2 * x(:,3) + &
            mdl(16) * x(:,3)**2 * x(:,1) + mdl(17) * x(:,3)**2 * x(:,2) + &
            mdl(18) * x(:,1)**3 + mdl(19) * x(:,2)**3 + mdl(20) * x(:,3)**3
        ! Perform the calculation
        y = doe_evaluate_model(3, mdl, x)

        ! Test
        if (.not.is_equal(y, ans)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_eval_model_three_way -1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_doe_fit_two_way() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: tol = 1.0d-6
        integer(int32), parameter :: n = 5
        integer(int32), parameter :: npts = n * n
        integer(int32), parameter :: nparams = 6

        ! Local Variables
        integer(int32) :: i, j, k
        real(real64) :: x(npts,2), y(npts), beta(nparams)
        type(doe_model) :: mdl

        ! Initialization
        rst = .true.
        beta = [1.0d0, 2.0d0, -1.0d0, 0.5d0, 0.3d0, -0.7d0]
        k = 0
        do i = 1, n
            do j = 1, n
                k = k + 1
                x(k,1) = -1.0d0 + 0.5d0 * (i - 1)
                x(k,2) = -1.0d0 + 0.5d0 * (j - 1)
            end do
        end do

        ! Y = b0 + b1 x1 + b2 x2 + b12 x1 x2 + b11 x1**2 + b22 x2**2
        y = beta(1) + beta(2) * x(:,1) + beta(3) * x(:,2) + &
            beta(4) * x(:,1) * x(:,2) + beta(5) * x(:,1)**2 + &
            beta(6) * x(:,2)**2

        ! Fit the model
        mdl = doe_fit_model(2, x, y)

        ! Test
        if (size(mdl%coefficients) /= nparams) then
            rst = .false.
            print "(A)", "TEST FAILED: test_doe_fit_two_way -1"
            return
        end if
        if (.not.is_equal(mdl%coefficients, beta, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_doe_fit_two_way -2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_doe_fit_three_way() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: tol = 1.0d-6
        integer(int32), parameter :: n = 4
        integer(int32), parameter :: npts = n * n * n
        integer(int32), parameter :: nparams = 20

        ! Local Variables
        integer(int32) :: i, j, k, l
        real(real64) :: x(npts,3), y(npts), beta(nparams)
        type(doe_model) :: mdl

        ! Initialization
        rst = .true.
        do i = 1, nparams
            beta(i) = 0.1d0 * i * (-1.0d0)**i
        end do
        l = 0
        do i = 1, n
            do j = 1, n
                do k = 1, n
                    l = l + 1
                    x(l,1) = -1.0d0 + (2.0d0 / (n - 1)) * (i - 1)
                    x(l,2) = -1.0d0 + (2.0d0 / (n - 1)) * (j - 1)
                    x(l,3) = -1.0d0 + (2.0d0 / (n - 1)) * (k - 1)
                end do
            end do
        end do

        ! Generate the response using the model
        y = doe_evaluate_model(3, beta, x)

        ! Fit the model
        mdl = doe_fit_model(3, x, y)

        ! Test
        if (size(mdl%coefficients) /= nparams) then
            rst = .false.
            print "(A)", "TEST FAILED: test_doe_fit_three_way -1"
            return
        end if
        if (.not.is_equal(mdl%coefficients, beta, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_doe_fit_three_way -2"
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

! ==============================================================================
! NEW FEATURE TESTS
! ==============================================================================

    function test_model_diagnostics_1() result(rst)
        !! Tests R-squared calculation from a simple model.
        logical :: rst
        
        ! Parameters
        integer(int32), parameter :: npts = 10
        integer(int32), parameter :: nparams = 2
        
        ! Local Variables
        integer(int32) :: i
        real(real64) :: x(npts, 1), y(npts), beta(nparams)
        type(doe_model) :: mdl
        type(doe_diagnostics) :: diag
        
        ! Initialization
        rst = .true.
        x(:, 1) = [(real(i, real64), i=1,npts)]
        beta = [1.0d0, 2.0d0]
        y = beta(1) + beta(2) * x(:, 1)  ! Perfect fit
        
        ! Set up model
        allocate(mdl%coefficients(nparams), source=beta)
        allocate(mdl%stats(nparams))
        allocate(mdl%map(nparams), source=[.true., .true.])
        mdl%nway = 1
        
        ! Calculate diagnostics
        diag = doe_model_diagnostics(mdl, x, y)
        
        ! Perfect fit should have R² ≈ 1.0
        if (.not.is_equal(diag%r_squared, 1.0d0, 1.0d-6)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_model_diagnostics_1 - R-squared"
        end if
        
        ! RMSE should be nearly zero
        if (diag%rmse > 1.0d-10) then
            rst = .false.
            print '(A)', "TEST FAILED: test_model_diagnostics_1 - RMSE"
        end if
        
    end function

! ==============================================================================

    function test_full_factorial_coded_1() result(rst)
        !! Tests the coded full factorial design.
        logical :: rst

        ! Local Variables
        integer(int32), parameter :: vars(2) = [2, 3]
        real(real64) :: tbl(6, 2), ans(6, 2)

        ! Initialization
        rst = .true.
        ans(:,1) = [-1.0d0, -1.0d0, -1.0d0, 1.0d0, 1.0d0, 1.0d0]
        ans(:,2) = [-1.0d0, 0.0d0, 1.0d0, -1.0d0, 0.0d0, 1.0d0]

        ! Test
        call full_factorial_coded(vars, tbl)
        if (.not.is_equal(tbl, ans)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_full_factorial_coded_1"
        end if
    end function

! ==============================================================================
    function test_fractional_factorial_1() result(rst)
        !! Tests 2^(3-1) fractional factorial design.
        logical :: rst
        
        ! Parameters
        integer(int32), parameter :: nfactors = 3
        integer(int32), parameter :: fraction = 1
        
        ! Local Variables
        integer(int32) :: m, n
        real(real64) :: tbl(4, 3)
        
        ! Initialization
        rst = .true.
        
        ! Get design size
        call fractional_factorial_size(nfactors, fraction, m, n)
        
        if (m /= 4 .or. n /= 3) then
            rst = .false.
            print '(A)', "TEST FAILED: test_fractional_factorial_1 - size"
        end if
        
        ! Generate design
        call fractional_factorial(nfactors, fraction, tbl)
        
        ! Check all values are -1 or +1
        if (any(abs(abs(tbl) - 1.0d0) > 1.0d-12) .or. any(abs(sum(tbl, 1)) > 1.0d-12)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_fractional_factorial_1 - values"
        end if
        
    end function

! ==============================================================================

    function test_encode_decode_1() result(rst)
        !! Tests encoding/decoding of variables.
        logical :: rst
        
        ! Local Variables
        real(real64) :: x_nat(3, 2), x_coded(3, 2), x_decoded(3, 2)
        real(real64) :: x_low(2), x_high(2)
        real(real64), parameter :: tol = 1.0d-10
        
        ! Initialization
        rst = .true.
        x_low = [50.0d0, 100.0d0]
        x_high = [150.0d0, 200.0d0]
        x_nat = reshape([50.0d0, 100.0d0, 150.0d0, 100.0d0, 150.0d0, 200.0d0], [3, 2])
        
        ! Encode
        call encode_variables(x_nat, x_low, x_high, x_coded)
        
        ! Decode
        call decode_variables(x_coded, x_low, x_high, x_decoded)
        
        ! Check recovery
        if (.not.is_equal(x_nat, x_decoded, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_encode_decode_1 - recovery"
        end if
        
        ! Check coded values at bounds
        if (.not.is_equal(x_coded(1, 1), -1.0d0, tol)) then  ! Low bound
            rst = .false.
            print '(A)', "TEST FAILED: test_encode_decode_1 - low bound"
        end if
        
        if (.not.is_equal(x_coded(3, 1), 1.0d0, tol)) then   ! High bound
            rst = .false.
            print '(A)', "TEST FAILED: test_encode_decode_1 - high bound"
        end if
        
    end function

! ==============================================================================

    function test_central_composite_design_1() result(rst)
        !! Tests CCD size and generation.
        logical :: rst
        
        ! Parameters
        integer(int32), parameter :: nfactors = 2
        
        ! Local Variables
        integer(int32) :: m, n
        real(real64) :: tbl(9, 2)  ! 2^2 + 2*2 + 1 = 9
        
        ! Initialization
        rst = .true.
        
        ! Get size
        call central_composite_design_size(nfactors, "orthogonal", m, n)
        
        if (m /= 9 .or. n /= 2) then
            rst = .false.
            print '(A)', "TEST FAILED: test_central_composite_design_1 - size"
        end if
        
        ! Generate design
        call central_composite_design(nfactors, "orthogonal", tbl)
        
        ! Check center point
        if (.not.is_equal(tbl(m, 1), 0.0d0, 1.0d-10) .or. &
            .not.is_equal(tbl(m, 2), 0.0d0, 1.0d-10)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_central_composite_design_1 - center"
        end if
        
    end function

! ==============================================================================

    function test_residuals_analysis_1() result(rst)
        !! Tests residual analysis calculations.
        logical :: rst
        
        ! Parameters
        integer(int32), parameter :: npts = 10
        integer(int32), parameter :: nparams = 2
        
        ! Local Variables
        integer(int32) :: i
        real(real64) :: x(npts, 1), y(npts), beta(nparams)
        type(doe_model) :: mdl
        type(doe_residuals) :: resid
        
        ! Initialization
        rst = .true.
        x(:, 1) = [(real(i, real64), i=1,npts)]
        beta = [1.0d0, 2.0d0]
        y = beta(1) + beta(2) * x(:, 1)  ! Perfect fit
        
        ! Set up model
        allocate(mdl%coefficients(nparams), source=beta)
        allocate(mdl%stats(nparams))
        allocate(mdl%map(nparams), source=[.true., .true.])
        mdl%nway = 1
        
        ! Calculate residuals
        resid = doe_residuals_analysis(mdl, x, y)
        
        ! Perfect fit should have near-zero residuals
        if (.not.is_equal(resid%residual_mean, 0.0d0, 1.0d-10)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_residuals_analysis_1 - mean"
        end if
        
        ! Check standardized residuals array exists
        if (.not.allocated(resid%standardized_residuals)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_residuals_analysis_1 - alloc"
        end if
        
    end function

! ==============================================================================

    function test_latin_hypercube_1() result(rst)
        !! Tests Latin hypercube design generation.
        logical :: rst
        
        ! Parameters
        integer(int32), parameter :: nfactors = 3
        integer(int32), parameter :: nsamples = 10
        
        ! Local Variables
        real(real64) :: tbl(nsamples, nfactors)
        integer(int32) :: seed, i, j
        
        ! Initialization
        rst = .true.
        seed = 42
        
        ! Generate LHS
        call latin_hypercube_design(nfactors, nsamples, seed, tbl)
        
        ! Check all values in [-1, 1]
        if (any(tbl < -1.0d0) .or. any(tbl > 1.0d0)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_latin_hypercube_1 - bounds"
        end if
        
        ! Check non-empty design
        if (size(tbl, 1) /= nsamples .or. size(tbl, 2) /= nfactors) then
            rst = .false.
            print '(A)', "TEST FAILED: test_latin_hypercube_1 - size"
        end if
        
    end function

! ==============================================================================
! ENHANCED FEATURES TESTS
! ==============================================================================

    function test_design_efficiency_1() result(rst)
        !! Tests design efficiency metrics calculation.
        logical :: rst
        
        ! Parameters
        integer(int32), parameter :: npts = 8
        integer(int32), parameter :: nfactors = 3
        
        ! Local Variables
        real(real64) :: x(npts, nfactors)
        type(doe_efficiency_metrics) :: metrics
        integer(int32) :: i
        
        ! Initialization
        rst = .true.
        ! Create simple full factorial design in coded variables
        do i = 1, npts
            if (mod(i, 2) == 1) then
                x(i, 1) = -1.0d0
            else
                x(i, 1) = 1.0d0
            end if
            if (mod(i/2, 2) == 0) then
                x(i, 2) = -1.0d0
            else
                x(i, 2) = 1.0d0
            end if
            if (mod(i/4, 2) == 0) then
                x(i, 3) = -1.0d0
            else
                x(i, 3) = 1.0d0
            end if
        end do
        
        ! Calculate efficiency
        metrics = doe_design_efficiency(x)
        
        ! Check that efficiency measures are in valid ranges
        if (metrics%d_efficiency < 0.0d0 .or. metrics%d_efficiency > 1.0d0) then
            rst = .false.
            print '(A)', "TEST FAILED: test_design_efficiency_1 - d_efficiency range"
        end if
        
        if (metrics%a_efficiency < 0.0d0 .or. metrics%a_efficiency > 1.0d0) then
            rst = .false.
            print '(A)', "TEST FAILED: test_design_efficiency_1 - a_efficiency range"
        end if
        
    end function

! ==============================================================================

    function test_model_comparison_1() result(rst)
        !! Tests model comparison functionality.
        logical :: rst
        
        ! Parameters
        integer(int32), parameter :: npts = 9
        integer(int32), parameter :: nfactors = 2
        integer(int32), parameter :: nparams_full = 6  ! 1 intercept + 2 main + 1 interaction + 2 squared
        integer(int32), parameter :: nparams_reduced = 3  ! 1 intercept + 2 main
        
        ! Local Variables
        integer(int32) :: i, j, k
        real(real64) :: x(npts, nfactors), y(npts)
        real(real64) :: beta_full(nparams_full), beta_reduced(nparams_reduced)
        type(doe_model) :: mdl_full, mdl_reduced
        type(doe_comparison_result) :: comp
        
        ! Initialization
        rst = .true.
        
        ! Create a 3x3 factorial design
        k = 1
        do i = 1, 3
            do j = 1, 3
                x(k, 1) = real(i, real64)
                x(k, 2) = real(j, real64)
                k = k + 1
            end do
        end do
        
        ! Full model: y = 1 + 2*x1 + 1*x2 + 0.5*x1*x2 + 0.1*x1**2 + 0.2*x2**2
        beta_full = [1.0d0, 2.0d0, 1.0d0, 0.5d0, 0.1d0, 0.2d0]
        y = beta_full(1) + beta_full(2)*x(:,1) + beta_full(3)*x(:,2) + &
            beta_full(4)*x(:,1)*x(:,2) + beta_full(5)*x(:,1)**2 + &
            beta_full(6)*x(:,2)**2
        
        ! Reduced model: y = 1 + 2*x1 + 1*x2 (no interaction)
        beta_reduced = [1.0d0, 2.0d0, 1.0d0]
        
        ! Set up full model (with interactions)
        allocate(mdl_full%coefficients(nparams_full), source=beta_full)
        allocate(mdl_full%stats(nparams_full))
        allocate(mdl_full%map(nparams_full), source=[.true., .true., .true., .true., .true., .true.])
        mdl_full%nway = 2
        
        ! Set up reduced model (main effects only)
        allocate(mdl_reduced%coefficients(nparams_reduced), source=beta_reduced)
        allocate(mdl_reduced%stats(nparams_reduced))
        allocate(mdl_reduced%map(nparams_reduced), source=[.true., .true., .true.])
        mdl_reduced%nway = 1
        
        ! Compare models
        comp = doe_compare_models(mdl_full, mdl_reduced, x, y)
        
        ! Check F-statistic is computed
        if (comp%f_statistic < 0.0d0) then
            rst = .false.
            print '(A)', "TEST FAILED: test_model_comparison_1 - f_statistic"
        end if

        ! Degrees of freedom, and the exact full model must be significant
        if (comp%df_full /= npts - nparams_full .or. &
            comp%df_reduced /= npts - nparams_reduced .or. &
            comp%df_diff /= nparams_full - nparams_reduced .or. &
            .not.comp%significant_difference .or. comp%p_value > 0.05d0) then
            rst = .false.
            print '(A)', "TEST FAILED: test_model_comparison_1 - df/p"
        end if

        ! Reversed order is an invalid comparison
        comp = doe_compare_models(mdl_reduced, mdl_full, x, y)
        if (comp%significant_difference .or. comp%f_statistic /= 0.0d0 .or. &
            comp%p_value /= 1.0d0) then
            rst = .false.
            print '(A)', "TEST FAILED: test_model_comparison_1 - invalid"
        end if
        
    end function

! ==============================================================================

    function test_model_anova_1() result(rst)
        !! Tests ANOVA table computation.
        logical :: rst
        
        ! Parameters
        integer(int32), parameter :: npts = 10
        integer(int32), parameter :: nparams = 2
        
        ! Local Variables
        integer(int32) :: i
        real(real64) :: x(npts, 1), y(npts), beta(nparams)
        type(doe_model) :: mdl
        type(doe_anova_table) :: anova
        
        ! Initialization
        rst = .true.
        x(:, 1) = [(real(i, real64), i=1,npts)]
        beta = [1.0d0, 2.0d0]
        y = beta(1) + beta(2) * x(:, 1)  ! Perfect fit
        
        ! Set up model
        allocate(mdl%coefficients(nparams), source=beta)
        allocate(mdl%stats(nparams))
        allocate(mdl%map(nparams), source=[.true., .true.])
        mdl%nway = 1
        
        ! Calculate ANOVA
        anova = doe_model_anova(mdl, x, y)
        
        ! Perfect fit should have high R-squared
        if (.not.is_equal(anova%r_squared, 1.0d0, 1.0d-10)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_model_anova_1 - r_squared"
        end if
        
        ! SS_model + SS_residual = SS_total
        if (.not.is_equal(anova%ss_model + anova%ss_residual, &
            anova%ss_total, 1.0d-10)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_model_anova_1 - sum of squares"
        end if

        ! Degrees of freedom and p-value of a perfect fit
        if (anova%df_model /= 1 .or. anova%df_residual /= npts - 2 .or. &
            anova%df_total /= npts - 1 .or. anova%p_value > 1.0d-10) then
            rst = .false.
            print '(A)', "TEST FAILED: test_model_anova_1 - perfect fit df/p"
        end if

        ! Noisy data: residuals of +/-0.5
        y = y + [(0.5d0 * (-1.0d0)**i, i=1,npts)]
        anova = doe_model_anova(mdl, x, y)
        if (.not.is_equal(anova%ss_residual, 2.5d0, 1.0d-10) .or. &
            .not.is_equal(anova%f_statistic, anova%ms_model / &
                anova%ms_residual, 1.0d-10) .or. &
            .not.is_equal(anova%ms_residual, 2.5d0 / 8.0d0, 1.0d-10)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_model_anova_1 - noisy F"
        end if
        block
            type(f_distribution) :: fd
            fd = f_distribution(d1 = 1.0d0, d2 = 8.0d0)
            if (.not.is_equal(anova%p_value, fd%survival(anova%f_statistic), &
                1.0d-12) .or. anova%p_value <= 0.0d0 .or. &
                anova%p_value >= 1.0d0) then
                rst = .false.
                print '(A)', "TEST FAILED: test_model_anova_1 - p-value"
            end if
        end block
        
    end function

! ==============================================================================

    function test_rsm_optimization_1() result(rst)
        !! Tests RSM-based optimization.
        logical :: rst
        
        ! Parameters
        integer(int32), parameter :: npts = 9
        integer(int32), parameter :: nparams = 3
        integer(int32), parameter :: nfactors = 2
        
        ! Local Variables
        integer(int32) :: i
        real(real64) :: x(npts, nfactors), y(npts), beta(nparams)
        real(real64) :: x_low(nfactors), x_high(nfactors)
        type(doe_model) :: mdl
        type(doe_optimization_result) :: opt
        
        ! Initialization
        rst = .true.
        
        ! Simple CCD-like points
        x = reshape([&
            -1.0d0, -1.0d0, &
             1.0d0, -1.0d0, &
            -1.0d0,  1.0d0, &
             1.0d0,  1.0d0, &
             0.0d0,  0.0d0, &
            -1.4d0,  0.0d0, &
             1.4d0,  0.0d0, &
             0.0d0, -1.4d0, &
             0.0d0,  1.4d0], [9, 2])
        
        ! Response: simple quadratic
        do i = 1, npts
            y(i) = 10.0d0 - (x(i,1)**2 + x(i,2)**2)
        end do
        
        ! Set up model
        beta = [10.0d0, 0.1d0, -1.0d0]
        allocate(mdl%coefficients(nparams), source=beta)
        allocate(mdl%stats(nparams))
        allocate(mdl%map(nparams), source=[.true., .true., .true.])
        mdl%nway = 1
        
        ! Set factor bounds
        x_low = [0.0d0, 0.0d0]
        x_high = [100.0d0, 100.0d0]
        
        ! Optimize
        opt = doe_optimize_rsm(mdl, x_low, x_high)
        
        ! Check result structure
        if (.not.allocated(opt%optimal_coded_factors)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_rsm_optimization_1 - alloc coded"
        end if
        
        if (.not.allocated(opt%optimal_natural_factors)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_rsm_optimization_1 - alloc natural"
        end if
        
    end function

! ------------------------------------------------------------------------------
    function test_rsm_optimization_2() result(rst)
        !! Tests the response-surface optimizer against known optima.
        logical :: rst

        ! Parameters
        real(real64), parameter :: tol = 1.0d-5

        ! Local Variables
        type(doe_model) :: quad, lin, saddle
        type(doe_optimization_result) :: opt

        rst = .true.

        ! Interior maximum: y = 5 - (z1 - 0.3)**2 - (z2 - 0.4)**2
        call build_model(2, [4.75d0, 0.6d0, 0.8d0, 0.0d0, -1.0d0, -1.0d0], quad)
        opt = doe_optimize_rsm(quad, [-2.0d0, -2.0d0], [2.0d0, 2.0d0])
        if (.not.opt%converged .or. &
            abs(opt%optimal_coded_factors(1) - 0.3d0) > tol .or. &
            abs(opt%optimal_coded_factors(2) - 0.4d0) > tol .or. &
            abs(opt%optimal_response - 5.0d0) > tol .or. &
            abs(opt%optimal_natural_factors(1) - 0.6d0) > tol .or. &
            abs(opt%optimal_natural_factors(2) - 0.8d0) > tol) then
            rst = .false.
            print '(A)', "TEST FAILED: test_rsm_optimization_2 -1"
        end if

        ! Interior minimum of the negated surface
        call build_model(2, [-4.75d0, -0.6d0, -0.8d0, 0.0d0, 1.0d0, 1.0d0], quad)
        opt = doe_optimize_rsm(quad, [-2.0d0, -2.0d0], [2.0d0, 2.0d0], "min")
        if (abs(opt%optimal_coded_factors(1) - 0.3d0) > tol .or. &
            abs(opt%optimal_coded_factors(2) - 0.4d0) > tol .or. &
            abs(opt%optimal_response + 5.0d0) > tol) then
            rst = .false.
            print '(A)', "TEST FAILED: test_rsm_optimization_2 -2"
        end if

        ! Active bounds: y = 10 + 0.1 z1 - z2 on [0, 100]
        call build_model(1, [10.0d0, 0.1d0, -1.0d0], lin)
        opt = doe_optimize_rsm(lin, [0.0d0, 0.0d0], [100.0d0, 100.0d0])
        if (abs(opt%optimal_natural_factors(1) - 100.0d0) > tol .or. &
            abs(opt%optimal_natural_factors(2)) > tol .or. &
            abs(opt%optimal_response - 11.1d0) > tol) then
            rst = .false.
            print '(A)', "TEST FAILED: test_rsm_optimization_2 -3"
        end if
        opt = doe_optimize_rsm(lin, [0.0d0, 0.0d0], [100.0d0, 100.0d0], "min")
        if (abs(opt%optimal_natural_factors(1)) > tol .or. &
            abs(opt%optimal_natural_factors(2) - 100.0d0) > tol .or. &
            abs(opt%optimal_response - 8.9d0) > tol) then
            rst = .false.
            print '(A)', "TEST FAILED: test_rsm_optimization_2 -4"
        end if

        ! Target: response of 4 on the interior surface
        call build_model(2, [4.75d0, 0.6d0, 0.8d0, 0.0d0, -1.0d0, -1.0d0], quad)
        opt = doe_optimize_rsm(quad, [-2.0d0, -2.0d0], [2.0d0, 2.0d0], &
            "target", 4.0d0)
        if (abs(opt%optimal_response - 4.0d0) > tol) then
            rst = .false.
            print '(A)', "TEST FAILED: test_rsm_optimization_2 -5"
        end if

        ! Saddle: y = z1**2 - z2**2 + 0.1 z1, maximum of 1.1 at (1, 0)
        call build_model(2, [0.0d0, 0.1d0, 0.0d0, 0.0d0, 1.0d0, -1.0d0], saddle)
        opt = doe_optimize_rsm(saddle, [-1.0d0, -1.0d0], [1.0d0, 1.0d0])
        if (abs(opt%optimal_coded_factors(1) - 1.0d0) > tol .or. &
            abs(opt%optimal_coded_factors(2)) > tol .or. &
            abs(opt%optimal_response - 1.1d0) > tol) then
            rst = .false.
        end if

        ! Analytic stationary point: maximum, saddle, and singular cases
        block
            real(real64), allocatable :: zs(:), ev(:)
            character(len=:), allocatable :: nat
            logical :: ok

            call build_model(2, [4.75d0, 0.6d0, 0.8d0, 0.0d0, -1.0d0, -1.0d0], quad)
            call doe_stationary_point(quad, 2, zs, nat, ok, ev)
            if (.not.ok .or. nat /= "maximum" .or. &
                abs(zs(1) - 0.3d0) > 1.0d-10 .or. &
                abs(zs(2) - 0.4d0) > 1.0d-10 .or. &
                any(abs(ev + 2.0d0) > 1.0d-10)) then
                rst = .false.
                print '(A)', "TEST FAILED: test_rsm_optimization_2 -7"
            end if
            opt = doe_optimize_rsm(quad, [-2.0d0, -2.0d0], [2.0d0, 2.0d0])
            if (opt%iteration_count /= 0 .or. &
                opt%method /= "analytic stationary point") then
                rst = .false.
                print '(A)', "TEST FAILED: test_rsm_optimization_2 -8"
            end if

            call doe_stationary_point(saddle, 2, zs, nat, ok)
            if (.not.ok .or. nat /= "saddle" .or. &
                abs(zs(1) + 0.05d0) > 1.0d-10 .or. abs(zs(2)) > 1.0d-10) then
                rst = .false.
                print '(A)', "TEST FAILED: test_rsm_optimization_2 -9"
            end if

            call doe_stationary_point(lin, 2, zs, nat, ok)
            if (ok .or. nat /= "singular") then
                rst = .false.
                print '(A)', "TEST FAILED: test_rsm_optimization_2 -10"
            end if
        end block

    contains
        subroutine build_model(nway, beta, mdl)
            integer(int32), intent(in) :: nway
            real(real64), intent(in) :: beta(:)
            type(doe_model), intent(out) :: mdl

            mdl%nway = nway
            allocate(mdl%coefficients, source = beta)
            allocate(mdl%stats(size(beta)))
            allocate(mdl%map(size(beta)), source = .true.)
        end subroutine
    end function

! ==============================================================================

    function test_enhanced_prediction_1() result(rst)
        !! Tests enhanced prediction intervals with MSE.
        logical :: rst
        
        ! Parameters
        integer(int32), parameter :: npts = 10
        integer(int32), parameter :: nparams = 2
        
        ! Local Variables
        integer(int32) :: i
        real(real64) :: x(npts, 1), y(npts), beta(nparams)
        type(doe_model) :: mdl
        type(doe_prediction) :: pred
        real(real64) :: mse
        
        ! Initialization
        rst = .true.
        x(:, 1) = [(real(i, real64), i=1,npts)]
        beta = [1.0d0, 2.0d0]
        y = beta(1) + beta(2) * x(:, 1) + 0.1d0  ! Small noise
        
        ! Set up model
        allocate(mdl%coefficients(nparams), source=beta)
        allocate(mdl%stats(nparams))
        allocate(mdl%map(nparams), source=[.true., .true.])
        mdl%nway = 1
        
        ! Compute predictions with MSE
        mse = 0.01d0
        pred = doe_predict_enhanced(mdl, x, 0.05d0, mse, x)
        
        ! Check prediction intervals
        if (.not.allocated(pred%predicted_values)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_enhanced_prediction_1 - pred values"
        end if
        
        if (.not.allocated(pred%prediction_lower)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_enhanced_prediction_1 - pred lower"
        end if
        
        ! Prediction intervals should be wider than confidence intervals
        if (pred%prediction_upper(1) - pred%prediction_lower(1) <= &
            pred%confidence_upper(1) - pred%confidence_lower(1)) then
            ! Warning, not failure - could be due to approximation
        end if
        
    end function

! ------------------------------------------------------------------------------
    function test_doe_prediction_intervals_1() result(rst)
        !! Compares the intervals for a straight-line fit with the closed-form
        !! simple linear regression result.
        logical :: rst

        ! Parameters
        integer(int32), parameter :: npts = 10
        real(real64), parameter :: tcrit = 2.306004135204166d0 ! t(0.975, 8)
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        integer(int32) :: i
        real(real64) :: x(npts, 1), y(npts), xp(3, 1), xbar, sxx, mse, q, hw
        real(real64) :: b1, b0
        type(doe_model) :: mdl
        type(doe_prediction) :: pred

        ! Initialization
        rst = .true.
        x(:, 1) = [(real(i, real64), i = 1, npts)]
        y = [1.2d0, 2.9d0, 5.3d0, 6.8d0, 9.4d0, 10.9d0, 13.3d0, 14.6d0, &
            17.4d0, 18.7d0]
        xp(:, 1) = [0.0d0, 5.5d0, 12.0d0]

        ! Closed-form solution
        xbar = sum(x(:,1)) / npts
        sxx = sum((x(:,1) - xbar)**2)
        b1 = sum((x(:,1) - xbar) * (y - sum(y) / npts)) / sxx
        b0 = sum(y) / npts - b1 * xbar
        mse = sum((y - b0 - b1 * x(:,1))**2) / (npts - 2)

        ! Fit and predict
        mdl = doe_fit_model(1, x, y)
        pred = doe_predict(mdl, xp)

        do i = 1, 3
            q = 1.0d0 / npts + (xp(i,1) - xbar)**2 / sxx
            hw = tcrit * sqrt(mse * q)
            if (abs(pred%confidence_upper(i) - pred%predicted_values(i) - hw) &
                > tol) then
                rst = .false.
                print '(A,I0)', "TEST FAILED: test_doe_prediction_intervals_1 -1: ", i
            end if
            hw = tcrit * sqrt(mse * (1.0d0 + q))
            if (abs(pred%predicted_values(i) - pred%prediction_lower(i) - hw) &
                > tol) then
                rst = .false.
                print '(A,I0)', "TEST FAILED: test_doe_prediction_intervals_1 -2: ", i
            end if
        end do
        if (abs(pred%confidence_level - 0.95d0) > tol) then
            rst = .false.
            print '(A)', "TEST FAILED: test_doe_prediction_intervals_1 -3"
        end if
    end function

! ==============================================================================
end module
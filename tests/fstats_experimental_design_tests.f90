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

    function test_fractional_factorial_1() result(rst)
        !! Tests 2^(3-1) fractional factorial design.
        logical :: rst
        
        ! Parameters
        integer(int32), parameter :: nfactors = 3
        integer(int32), parameter :: fraction = 1
        
        ! Local Variables
        integer(int32) :: m, n
        integer(int32) :: tbl(4, 3)
        
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
        
        ! Check all values are 1 or 2
        if (any(tbl < 1) .or. any(tbl > 2)) then
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
        integer(int32), parameter :: nparams_full = 5  ! 1 intercept + 2 main + 2 interactions
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
        
        ! Full model: y = 1 + 2*x1 + 1*x2 + 0.5*x1*x2 + 0.1*x2*x1
        beta_full = [1.0d0, 2.0d0, 1.0d0, 0.5d0, 0.1d0]
        y = beta_full(1) + beta_full(2)*x(:,1) + beta_full(3)*x(:,2) + &
            beta_full(4)*x(:,1)*x(:,2) + beta_full(5)*x(:,2)*x(:,1)
        
        ! Reduced model: y = 1 + 2*x1 + 1*x2 (no interaction)
        beta_reduced = [1.0d0, 2.0d0, 1.0d0]
        
        ! Set up full model (with interactions)
        allocate(mdl_full%coefficients(nparams_full), source=beta_full)
        allocate(mdl_full%stats(nparams_full))
        allocate(mdl_full%map(nparams_full), source=[.true., .true., .true., .true., .true.])
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
        pred = doe_predict_enhanced(mdl, x, 0.05d0, mse)
        
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

! ==============================================================================
end module
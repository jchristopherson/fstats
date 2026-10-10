module fstats_glm_tests
    use iso_fortran_env
    use ieee_arithmetic
    use fstats, only : irls, glm_options, glm_convergence_info, &
        glm_coefficient_statistics, robust_weight_function, &
        FS_GLM_BINOMIAL, FS_GLM_POISSON, FS_GLM_GAMMA
    use fortran_test_helper
    implicit none
    private
    public :: test_glm_irls, test_glm_invalid_input
    character(32) :: callback_case

contains
    function test_glm_irls() result(rst)
        logical :: rst
        real(real64), dimension(4,1) :: x
        real(real64), dimension(4,2) :: design
        real(real64), dimension(4) :: y, mu
        real(real64), dimension(1) :: beta, saved
        real(real64), dimension(2) :: coeff, expected
        type(glm_options) :: opts, defaults, copied
        type(glm_convergence_info) :: info
        procedure(robust_weight_function), pointer :: weights
        integer(int32) :: family, scale_case
        real(real64) :: magnitude
        logical :: halt_overflow, halt_invalid, halt_divide

        call ieee_get_halting_mode(ieee_overflow, halt_overflow)
        call ieee_get_halting_mode(ieee_invalid, halt_invalid)
        call ieee_get_halting_mode(ieee_divide_by_zero, halt_divide)
        call ieee_set_halting_mode(ieee_overflow, .true.)
        call ieee_set_halting_mode(ieee_invalid, .true.)
        call ieee_set_halting_mode(ieee_divide_by_zero, .true.)
        rst = .true.
        call defaults%set_to_default()
        copied = opts
        if (copied%max_iteration_count /= defaults%max_iteration_count .or. &
            copied%tolerance /= defaults%tolerance .or. &
            (copied%use_robust_weighting .neqv. defaults%use_robust_weighting)) then
            rst = .false.
            print '(A)', 'TEST FAILED: GLM initialized defaults'
        end if
        if (copied%robust_weighting_constant /= defaults%robust_weighting_constant) then
            rst = .false.
            print '(A)', 'TEST FAILED: GLM initialized tuning constant'
        end if
        opts%use_robust_weighting = .false.
        x = 1.0d0
        do family = FS_GLM_BINOMIAL, FS_GLM_GAMMA
            select case (family)
            case (FS_GLM_BINOMIAL)
                y = [0.0d0, 0.0d0, 0.0d0, 1.0d0]
                expected(1) = log(1.0d0 / 3.0d0)
                beta = 0.0d0
            case (FS_GLM_POISSON)
                y = [0.0d0, 1.0d0, 2.0d0, 5.0d0]
                expected(1) = log(2.0d0)
                beta = 0.0d0
            case (FS_GLM_GAMMA)
                y = [1.0d0, 2.0d0, 3.0d0, 6.0d0]
                expected(1) = 1.0d0 / 3.0d0
                beta = 2.0d0
            end select
            call irls(family, x, y, beta, options = opts, info = info)
            if (.not.info%converged .or. .not.assert(beta(1), expected(1), 1.0d-8) .or. &
                info%iteration_count < 1 .or. info%residual_value >= opts%tolerance) then
                rst = .false.
                print '(A,I0)', 'TEST FAILED: GLM intercept family ', family
            end if
        end do

        design(:,1) = 1.0d0
        design(:,2) = [-1.0d0, -1.0d0, 1.0d0, 1.0d0]
        do family = FS_GLM_BINOMIAL, FS_GLM_GAMMA
            coeff = 0.0d0
            select case (family)
            case (FS_GLM_BINOMIAL)
                y = [0.1d0, 0.3d0, 0.7d0, 0.9d0]
                expected = [0.0d0, log(4.0d0)]
            case (FS_GLM_POISSON)
                y = [1.0d0, 3.0d0, 4.0d0, 8.0d0]
                expected = [0.5d0 * log(12.0d0), 0.5d0 * log(3.0d0)]
            case (FS_GLM_GAMMA)
                y = [1.0d0, 3.0d0, 4.0d0, 8.0d0]
                expected = [1.0d0 / 3.0d0, -1.0d0 / 6.0d0]
                coeff(1) = 0.2d0
            end select
            call irls(family, design, y, coeff, options = opts, info = info)
            if (.not.info%converged .or. .not.assert(coeff, expected, 1.0d-8)) then
                rst = .false.
                print '(A,I0)', 'TEST FAILED: GLM two-column family ', family
            end if
            mu = matmul(design, coeff)
            select case (family)
            case (FS_GLM_BINOMIAL)
                mu = 1.0d0 / (1.0d0 + exp(-mu))
            case (FS_GLM_POISSON)
                mu = exp(mu)
            case (FS_GLM_GAMMA)
                mu = 1.0d0 / mu
            end select
            if (maxval(abs(matmul(transpose(design), y - mu))) > 1.0d-8) then
                rst = .false.
                print '(A,I0)', 'TEST FAILED: GLM score equations family ', family
            end if
        end do

        x = 1.0d0
        do scale_case = 1, 4
            magnitude = 1.0d200
            if (scale_case == 2) magnitude = 1.0d-200
            if (scale_case == 3) magnitude = 1.0d308
            if (scale_case == 4) magnitude = 1.0d-308
            y = magnitude
            beta = log(magnitude)
            opts%tolerance = 1.0d-8
            call irls(FS_GLM_POISSON, x, y, beta, options = opts, info = info)
            if (.not.info%converged .or. abs(beta(1) - log(magnitude)) > 1.0d-8) then
                rst = .false.
                print '(A)', 'TEST FAILED: GLM extreme Poisson mean'
            end if
            beta = 1.0d0 / magnitude
            opts%tolerance = abs(beta(1)) * 1.0d-10
            call irls(FS_GLM_GAMMA, x, y, beta, options = opts, info = info)
            if (.not.info%converged .or. abs(beta(1) * magnitude - 1.0d0) > 1.0d-8) then
                rst = .false.
                print '(A)', 'TEST FAILED: GLM extreme Gamma mean'
            end if
            x = magnitude
            y = [0.0d0, 1.0d0, 2.0d0, 5.0d0]
            beta = 0.0d0
            opts%tolerance = 1.0d-10 / magnitude
            call irls(FS_GLM_POISSON, x, y, beta, options = opts, info = info)
            if (.not.info%converged .or. abs(beta(1) * magnitude - log(2.0d0)) > 1.0d-8) then
                rst = .false.
                print '(A)', 'TEST FAILED: GLM extreme design scale'
            end if
            x = 1.0d0
        end do
        design(:,1) = 1.0d0
        design(:,2) = 1.0d0 + 1.0d-7 * [-1.0d0, -1.0d0, 1.0d0, 1.0d0]
        y = exp([0.8d0, 0.8d0, 1.2d0, 1.2d0])
        coeff = 0.0d0
        opts%tolerance = 1.0d-2
        call irls(FS_GLM_POISSON, design, y, coeff, options = opts, info = info)
        if (.not.info%converged .or. &
            maxval(abs(matmul(design, coeff) - log(y))) > 1.0d-8) then
            rst = .false.
            print '(A)', 'TEST FAILED: GLM nearly collinear design'
        end if
        design(:,1) = 1.0d200
        design(:,2) = 1.0d-200 * [-1.0d0, -1.0d0, 1.0d0, 1.0d0]
        coeff = 0.0d0
        opts%tolerance = 1.0d190
        call irls(FS_GLM_POISSON, design, y, coeff, options = opts, info = info)
        if (.not.info%converged .or. &
            maxval(abs(matmul(design, coeff) - log(y))) > 1.0d-8) then
            rst = .false.
            print '(A)', 'TEST FAILED: GLM independently scaled columns'
        end if
        opts%tolerance = 1.0d-8
        beta = -500.0d0
        y = exp(-500.0d0)
        call irls(FS_GLM_BINOMIAL, x, y, beta, options = opts, info = info)
        if (.not.info%converged .or. beta(1) /= -500.0d0) then
            rst = .false.
            print '(A)', 'TEST FAILED: GLM extreme binomial predictor'
        end if

        y = 2.0d0
        beta = -10.0d0
        call irls(FS_GLM_POISSON, x, y, beta, options = opts, info = info)
        if (.not.info%converged .or. abs(beta(1) - log(2.0d0)) > 1.0d-8) then
            rst = .false.
            print '(A)', 'TEST FAILED: GLM Poisson domain-preserving steps'
        end if

        beta = 0.0d0
        y = 4.0d0
        opts%max_iteration_count = 1
        call irls(FS_GLM_POISSON, x, y, beta, options = opts, info = info)
        if (info%converged .or. info%iteration_count /= 1 .or. &
            .not.assert(beta(1), 1.5d0, 1.0d-12) .or. &
            .not.assert(info%residual_value, 1.5d0, 1.0d-12)) then
            rst = .false.
            print '(A)', 'TEST FAILED: GLM iteration limit'
        end if
        opts%max_iteration_count = 100

        opts%use_robust_weighting = .true.
        callback_case = 'selected'
        weights => test_weights
        y = [1.0d0, 9.0d0, 3.0d0, 9.0d0]
        beta = 0.0d0
        call irls(FS_GLM_POISSON, x, y, beta, weights, opts, info)
        if (.not.info%converged .or. abs(beta(1) - log(2.0d0)) > 1.0d-8) then
            rst = .false.
            print '(A)', 'TEST FAILED: GLM custom zero/nonzero weights'
        end if
        y = [4.0d0, huge(1.0d0), 4.0d0, huge(1.0d0)]
        beta = 0.0d0
        opts%max_iteration_count = 1
        callback_case = 'selected'
        call irls(FS_GLM_POISSON, x, y, beta, weights, opts, info)
        if (info%converged .or. abs(beta(1) - 1.5d0) > 1.0d-12) then
            rst = .false.
            print '(A)', 'TEST FAILED: GLM zero-weight deviance scaling'
        end if
        opts%max_iteration_count = 100
        callback_case = 'huge'
        y = [0.0d0, 1.0d0, 2.0d0, 5.0d0]
        beta = 0.0d0
        call irls(FS_GLM_POISSON, x, y, beta, weights, opts, info)
        if (.not.info%converged .or. abs(beta(1) - log(2.0d0)) > 1.0d-8) then
            rst = .false.
            print '(A)', 'TEST FAILED: GLM huge robust weights'
        end if
        nullify(weights)
        beta = 0.0d0
        y = [0.0d0, 1.0d0, 1.0d0, 1.0d0]
        call irls(FS_GLM_BINOMIAL, x, y, beta, weights = weights, info = info)
        saved = beta
        beta = 0.0d0
        call irls(FS_GLM_BINOMIAL, x, y, beta, info = info)
        if (.not.info%converged .or. .not.assert(beta, saved, 1.0d-12)) then
            rst = .false.
            print '(A)', 'TEST FAILED: GLM absent/null callback defaults'
        end if
        callback_case = 'unallocated'
        weights => test_weights
        opts%use_robust_weighting = .false.
        opts%robust_weighting_constant = ieee_value(0.0d0, ieee_quiet_nan)
        beta = 0.0d0
        call irls(FS_GLM_BINOMIAL, x, y, beta, weights, opts, info)
        if (.not.info%converged .or. abs(beta(1) - log(3.0d0)) > 1.0d-8) then
            rst = .false.
            print '(A)', 'TEST FAILED: GLM disabled robust callback'
        end if
        if (.not.test_glm_inference()) rst = .false.
        if (rst) print '(A)', 'GLM numerical, fit, and inference tests passed.'
        call ieee_set_halting_mode(ieee_overflow, halt_overflow)
        call ieee_set_halting_mode(ieee_invalid, halt_invalid)
        call ieee_set_halting_mode(ieee_divide_by_zero, halt_divide)
    end function

    function test_glm_inference() result(rst)
        logical :: rst
        type(glm_options) :: opts
        type(glm_convergence_info) :: info
        type(glm_coefficient_statistics), dimension(1) :: stats
        type(glm_coefficient_statistics), dimension(2) :: two_stats
        real(real64), dimension(4,1) :: x
        real(real64), dimension(4,2) :: design
        real(real64), dimension(4) :: y
        real(real64), dimension(1) :: beta
        real(real64), dimension(2) :: coefficients
        real(real64), dimension(1,1) :: covariance, covariance_only
        real(real64), dimension(2,2) :: two_covariance, expected
        real(real64) :: variance, critical, wald, weight, derivative, cutoff
        integer(int32) :: family
        procedure(robust_weight_function), pointer :: weights
        rst = .true.
        opts%use_robust_weighting = .false.
        x = 1.0d0
        do family = FS_GLM_BINOMIAL, FS_GLM_GAMMA
            beta = 0.0d0
            critical = 1.959963984540054d0
            select case (family)
            case (FS_GLM_BINOMIAL)
                y = [0.0d0, 0.0d0, 0.0d0, 1.0d0]
                variance = 4.0d0 / 3.0d0
            case (FS_GLM_POISSON)
                y = [0.0d0, 1.0d0, 2.0d0, 5.0d0]
                variance = 0.125d0
            case (FS_GLM_GAMMA)
                y = [1.0d0, 2.0d0, 3.0d0, 6.0d0]
                beta = 0.2d0
                variance = (14.0d0 / 27.0d0) / 36.0d0
                critical = 3.182446305284263d0
            end select
            call irls(family, x, y, beta, options = opts, info = info, &
                stats = stats, cov = covariance)
            wald = beta(1) / sqrt(variance)
            if (.not.info%converged .or. .not.assert(covariance(1,1), variance, 1.0d-8) .or. &
                .not.assert(stats(1)%standard_error, sqrt(variance), 1.0d-8) .or. &
                .not.assert(stats(1)%wald_statistic, wald, 1.0d-8) .or. &
                .not.assert(stats(1)%confidence_interval_lower, beta(1) - critical * sqrt(variance), 1.0d-7) .or. &
                .not.assert(stats(1)%confidence_interval_upper, beta(1) + critical * sqrt(variance), 1.0d-7)) then
                rst = .false.
                print '(A,I0)', 'TEST FAILED: GLM analytic inference family ', family
            end if
            if (family /= FS_GLM_GAMMA) then
                if (.not.assert(stats(1)%p_value, erfc(abs(wald) / sqrt(2.0d0)), 1.0d-8)) then
                    rst = .false.
                    print '(A)', 'TEST FAILED: GLM normal Wald p-value'
                end if
            else
                if (abs(stats(1)%p_value - (1.0d0 - (2.0d0 / acos(-1.0d0)) * &
                    (atan(abs(wald) / sqrt(3.0d0)) + sqrt(3.0d0) * abs(wald) / &
                    (wald**2 + 3.0d0)))) > 1.0d-8) then
                    rst = .false.
                    print '(A)', 'TEST FAILED: GLM Gamma t p-value'
                end if
            end if
            call irls(family, x, y, beta, options = opts, cov = covariance_only)
            if (.not.assert(covariance_only, covariance, 1.0d-8)) then
                rst = .false.
                print '(A)', 'TEST FAILED: GLM covariance-only output'
            end if
        end do
        design(:,1) = 10.0d0
        design(:,2) = [0.0d0, 0.0d0, 0.1d0, 0.1d0]
        y = [1.0d0, 3.0d0, 4.0d0, 8.0d0]
        coefficients = 0.0d0
        expected(1,:) = [0.0025d0, -0.25d0]
        expected(2,:) = [-0.25d0, 100.0d0 / 3.0d0]
        call irls(FS_GLM_POISSON, design, y, coefficients, options = opts, &
            stats = two_stats, cov = two_covariance)
        if (.not.assert(two_covariance, expected, 1.0d-7) .or. &
            .not.assert(two_stats(1)%standard_error**2, two_covariance(1,1), 1.0d-10) .or. &
            .not.assert(two_stats(2)%standard_error**2, two_covariance(2,2), 1.0d-7)) then
            rst = .false.
            print '(A)', 'TEST FAILED: GLM scaled/pivoted full covariance'
        end if
        weights => test_weights
        callback_case = 'huge'
        opts%use_robust_weighting = .true.
        expected(2,2) = 325.0d0 / 9.0d0
        call irls(FS_GLM_POISSON, design, y, coefficients, weights = weights, options = opts, &
            cov = two_covariance)
        if (.not.assert(two_covariance, expected, 1.0d-6)) then
            rst = .false.
            print '(A)', 'TEST FAILED: GLM custom-weight full sandwich covariance'
        end if
        opts%use_robust_weighting = .false.
        opts%max_iteration_count = 100
        y = 1.0d308
        beta = log(1.0d308)
        call irls(FS_GLM_POISSON, x, y, beta, options = opts, stats = stats, cov = covariance)
        if (abs(covariance(1,1) * 1.0d308 - 0.25d0) > 1.0d-10 .or. &
            abs(stats(1)%standard_error * 1.0d154 - 0.5d0) > 1.0d-10) then
            rst = .false.
            print '(A)', 'TEST FAILED: GLM subnormal covariance'
        end if
        y = 1.0d-308
        beta = log(1.0d-308)
        call irls(FS_GLM_POISSON, x, y, beta, options = opts, stats = stats, cov = covariance)
        if (abs(covariance(1,1) * 1.0d-308 - 0.25d0) > 1.0d-10 .or. &
            abs(stats(1)%standard_error * 1.0d-154 - 0.5d0) > 1.0d-10) then
            rst = .false.
            print '(A)', 'TEST FAILED: GLM near-limit covariance'
        end if
        y = 1.0d0
        beta = 1.0d0
        call irls(FS_GLM_GAMMA, x, y, beta, options = opts, stats = stats, cov = covariance)
        if (covariance(1,1) /= 0.0d0 .or. stats(1)%standard_error /= 0.0d0 .or. &
            stats(1)%p_value /= 0.0d0 .or. .not.ieee_is_finite(stats(1)%confidence_interval_lower)) then
            rst = .false.
            print '(A)', 'TEST FAILED: GLM zero Pearson dispersion'
        end if
        beta = 0.0d0
        y = [0.0d0, 0.0d0, 1.0d0, 1.0d0]
        opts%use_robust_weighting = .true.
        cutoff = opts%robust_weighting_constant
        weight = (1.0d0 - (0.5d0 / cutoff)**2)**2
        derivative = -4.0d0 * 0.5d0 / cutoff**2 * (1.0d0 - (0.5d0 / cutoff)**2)
        variance = (4.0d0 / 3.0d0) * weight**2 / (weight + 0.5d0 * derivative)**2
        call irls(FS_GLM_BINOMIAL, x, y, beta, options = opts, stats = stats, cov = covariance)
        if (.not.assert(covariance(1,1), variance, 1.0d-7) .or. &
            .not.assert(stats(1)%p_value, 1.0d0, 1.0d-10)) then
            rst = .false.
            print '(A)', 'TEST FAILED: GLM differentiated Tukey sandwich'
        end if
        opts%use_robust_weighting = .false.
        opts%max_iteration_count = 1
        beta = 0.0d0
        y = 4.0d0
        call irls(FS_GLM_POISSON, x, y, beta, options = opts, stats = stats, cov = covariance)
        if (.not.all(ieee_is_nan(covariance)) .or. .not.ieee_is_nan(stats(1)%standard_error) .or. &
            .not.ieee_is_nan(stats(1)%wald_statistic) .or. .not.ieee_is_nan(stats(1)%p_value) .or. &
            .not.ieee_is_nan(stats(1)%confidence_interval_lower) .or. &
            .not.ieee_is_nan(stats(1)%confidence_interval_upper)) then
            rst = .false.
            print '(A)', 'TEST FAILED: GLM nonconvergence inference NaNs'
        end if
    end function

    subroutine test_glm_invalid_input(which)
        character(*), intent(in) :: which
        real(real64), dimension(4,2) :: x
        real(real64), dimension(4) :: y
        real(real64), dimension(2) :: beta
        real(real64), dimension(2,2) :: covariance
        type(glm_coefficient_statistics), dimension(2) :: stats
        type(glm_options) :: opts
        procedure(robust_weight_function), pointer :: weights
        integer(int32) :: family
        call ieee_set_halting_mode(ieee_overflow, .true.)
        call ieee_set_halting_mode(ieee_invalid, .true.)
        call ieee_set_halting_mode(ieee_divide_by_zero, .true.)
        x(:,1) = 1.0d0
        x(:,2) = [-1.0d0, -1.0d0, 1.0d0, 1.0d0]
        y = 0.5d0
        beta = 0.0d0
        family = FS_GLM_BINOMIAL
        opts%use_robust_weighting = .false.
        nullify(weights)
        select case (which)
        case ('stats-size')
            call irls(family, x, y, beta, options = opts, stats = stats(:1))
            return
        case ('covariance-size')
            call irls(family, x, y, beta, options = opts, cov = covariance(:1,:))
            return
        case ('alpha-zero')
            call irls(family, x, y, beta, alpha = 0.0d0)
            return
        case ('alpha-one')
            call irls(family, x, y, beta, alpha = 1.0d0)
            return
        case ('alpha-nan')
            call irls(family, x, y, beta, alpha = ieee_value(0.0d0, ieee_quiet_nan))
            return
        case ('gamma-inference-dof')
            call irls(FS_GLM_GAMMA, x(:2,:), y(:2), beta, options = opts, cov = covariance)
            return
        case ('robust-inference-dof')
            call irls(family, x(:2,:), y(:2), beta, cov = covariance)
            return
        case ('covariance-overflow')
            family = FS_GLM_POISSON
            x(:,1) = 1.0d-200
            y = 1.0d0
            call irls(family, x(:,:1), y, beta(:1), options = opts, cov = covariance(:1,:1))
            return
        case ('rows')
            call irls(family, x, y(:3), beta)
            return
        case ('columns')
            call irls(family, x, y, beta(:1))
            return
        case ('empty')
            call irls(family, x(:0,:), y(:0), beta)
            return
        case ('no-columns')
            call irls(family, x(:,:0), y, beta(:0))
            return
        case ('underdetermined')
            call irls(family, x(:1,:), y(:1), beta)
            return
        case ('family')
            family = 99
        case ('nan-x')
            x(1,1) = ieee_value(0.0d0, ieee_quiet_nan)
        case ('inf-y')
            y(1) = ieee_value(0.0d0, ieee_positive_inf)
        case ('nan-beta')
            beta(1) = ieee_value(0.0d0, ieee_quiet_nan)
        case ('binomial-range')
            y(1) = 1.1d0
        case ('poisson-range')
            family = FS_GLM_POISSON
            y(1) = -1.0d0
        case ('gamma-range')
            family = FS_GLM_GAMMA
            y(1) = 0.0d0
        case ('iterations')
            opts%max_iteration_count = 0
        case ('tolerance')
            opts%tolerance = 0.0d0
        case ('nan-tolerance')
            opts%tolerance = ieee_value(0.0d0, ieee_quiet_nan)
        case ('cutoff')
            opts%use_robust_weighting = .true.
            opts%robust_weighting_constant = -1.0d0
        case ('inf-cutoff')
            opts%use_robust_weighting = .true.
            opts%robust_weighting_constant = ieee_value(0.0d0, ieee_positive_inf)
        case ('singular')
            x(:,2) = x(:,1)
        case ('zero-column')
            x(:,2) = 0.0d0
        case ('gamma-domain')
            family = FS_GLM_GAMMA
        case ('poisson-overflow')
            family = FS_GLM_POISSON
            beta(1) = 1000.0d0
        case ('poisson-underflow')
            family = FS_GLM_POISSON
            beta(1) = -1000.0d0
        case ('binomial-boundary')
            beta(1) = 1000.0d0
        case ('predictor-overflow')
            x(:,1) = huge(1.0d0)
            beta(1) = 2.0d0
        case ('sum-overflow')
            x = 1.0d0
            beta = huge(1.0d0)
        case ('negative', 'nan', 'infinite', 'size', 'unallocated', 'zero')
            opts%use_robust_weighting = .true.
            callback_case = which
            weights => test_weights
        case default
            stop 0
        end select
        call irls(family, x, y, beta, weights, opts)
    end subroutine

    function test_weights(r, c) result(w)
        real(real64), intent(in), dimension(:) :: r
        real(real64), intent(in) :: c
        real(real64), allocatable, dimension(:) :: w
        if (callback_case == 'unallocated') return
        allocate(w(size(r)))
        w = 1.0d0
        select case (callback_case)
        case ('selected')
            w = [1.0d0, 0.0d0, 1.0d0, 0.0d0]
        case ('huge')
            w = huge(c)
        case ('negative')
            w(1) = -1.0d0
        case ('nan')
            w(1) = ieee_value(0.0d0, ieee_quiet_nan)
        case ('infinite')
            w(1) = ieee_value(0.0d0, ieee_positive_inf)
        case ('size')
            w = [1.0d0]
        case ('zero')
            w = 0.0d0
        end select
    end function
end module
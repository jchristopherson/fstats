module fstats_glm_irls
    !! Numerically scaled IRLS for logit-binomial, log-Poisson, and
    !! reciprocal-Gamma generalized linear models.
    use iso_fortran_env
    use ieee_arithmetic
    use fstats_glm_types
    use fstats_errors
    use fstats_robust_statistics
    use fstats_regression, only : iteration_update
    use fstats_normal_distribution, only : normal_distribution
    use fstats_t_distribution, only : t_distribution
    use fstats_hypothesis, only : confidence_interval
    use linalg, only : qr_factor, mult_qr, solve_triangular_system
    implicit none
    private
    public :: irls
    public :: glm_iteration_update
    public :: glm_binomial_link_inv
    public :: glm_poisson_link_inv
    public :: glm_gamma_link_inv

    abstract interface
        subroutine glm_iteration_update(iter, resid, beta)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: iter
            real(real64), intent(in), dimension(:) :: resid
            real(real64), intent(in), dimension(:) :: beta
        end subroutine
    end interface

contains
! ------------------------------------------------------------------------------
    subroutine irls(method, x, y, beta, weights, options, info, status, resid, &
        stats, alpha, cov)
        !! Fit a GLM from initial coefficients. The working response and
        !! weights are z = eta + (y - mu) * g'(mu) and
        !! W = 1 / (g'(mu)**2 * V(mu)). The implementation solves the equivalent
        !! increment problem sqrt(W) X delta = sqrt(W) (z - eta), using
        !! column-scaled, pivoted QR rather than normal equations. Weights and
        !! right-hand sides are normalized in log space to avoid overflow.
        !!
        !! No intercept, offset, binomial trial counts, or prior observation
        !! weights are added. A common dispersion factor cancels from updates.
        !! Robust weights act on raw residuals without scale estimation; set
        !! options%use_robust_weighting false for ordinary likelihood fitting.
        !!
        !! Initial predictors must be in the family's domain. Proposed updates
        !! are halved up to 64 times to preserve that domain and avoid increasing
        !! deviance (with current robust weights held fixed). Only an unhalved
        !! step can declare convergence. There is no
        !! guarantee of convergence, especially for separated binomial data.
        !! The iteration limit returns the last valid coefficients; supply info
        !! to distinguish convergence from an exhausted iteration budget.
        !!
        !! Size mismatches stop with FS_MATRIX_SIZE_ERROR. Empty or nonfinite
        !! inputs, invalid family/options/responses/callbacks, unrepresentable
        !! arithmetic, and invalid mean domains stop with FS_INVALID_INPUT_ERROR.
        !! Too few observations or numerical rank deficiency of the weighted,
        !! scaled design stop with FS_UNDERDEFINED_PROBLEM_ERROR. Missing data
        !! must be removed or imputed by the caller.
        integer(int32), intent(in) :: method
            !! FS_GLM_BINOMIAL (logit), FS_GLM_POISSON (log), or
            !! FS_GLM_GAMMA (reciprocal).
        real(real64), intent(in), dimension(:,:) :: x
            !! Finite N-by-P design, N >= P > 0. Weighted columns must be
            !! numerically independent at epsilon * max(N,P) relative scale.
        real(real64), intent(in), dimension(:) :: y
            !! N finite responses: [0,1] for binomial, >= 0 for Poisson,
            !! > 0 for Gamma. Fractional responses are allowed; binomial
            !! proportions are treated as unit-trial quasi-likelihood data.
        real(real64), intent(inout), dimension(:) :: beta
            !! P finite starting coefficients; receives the last valid iterate.
            !! Zero is valid for binomial/Poisson, but Gamma requires X beta > 0.
        procedure(robust_weight_function), pointer, intent(in), optional :: weights
            !! Robust callback used only when enabled. Absent/disassociated
            !! pointers select Tukey's biweight on raw response residuals.
        type(glm_options), intent(in), optional :: options
            !! Solver controls; omission uses the documented component defaults.
        type(glm_convergence_info), intent(out), optional :: info
            !! Convergence flag, completed update count, and last maximum
            !! absolute coefficient change (not the response residual).
        procedure(glm_iteration_update), intent(in), pointer, optional :: status
            !! An optional pointer to a routine that can be used to extract
            !! iteration information.
        real(real64), intent(out), optional, dimension(:) :: resid
            !! An optional P-element array where the residual will be written.
        type(glm_coefficient_statistics), intent(out), optional, dimension(:) :: stats
            !! P coefficient statistics, computed only after convergence;
            !! otherwise all fields are NaN. Ordinary binomial/Poisson fits use
            !! unit-dispersion information covariance and normal Wald inference.
            !! Ordinary Gamma fits use Pearson dispersion and approximate t
            !! inference with N-P degrees of freedom. Robust fits use an HC1
            !! empirical sandwich and normal Wald inference, treating rows as
            !! independent. Robust callbacks must be deterministic, smooth near
            !! the solution, and observation-wise for this inference to apply.
            !! Gamma and robust inference require N > P. Fractional responses
            !! do not automatically enable overdispersion estimation.
        real(real64), intent(in), optional :: alpha
            !! Finite significance level in (0,1), default 0.05. Intervals have
            !! confidence level 1-alpha. Validated even when stats is omitted.
        real(real64), intent(out), optional, dimension(:,:) :: cov
            !! P-by-P coefficient covariance in original beta order and units.
            !! Uses the same final-model information/Pearson or HC1 sandwich
            !! calculation as stats. May be requested independently of stats.
            !! All entries are NaN if fitting does not converge. Gamma and
            !! robust covariance require N > P; unrepresentable entries stop
            !! with FS_INVALID_INPUT_ERROR. Not a response covariance matrix.

        ! Local Variables
        procedure(inverse_link_function), pointer :: link_inv
        type(glm_convergence_info) :: fit_info
        type(glm_options) :: opts
        real(real64) :: significance, nan

        ! Input Validation
        if (size(x, 1) /= size(y) .or. size(x, 2) /= size(beta)) &
            error stop FS_MATRIX_SIZE_ERROR
        if (size(y) == 0 .or. size(beta) == 0) &
            error stop FS_INVALID_INPUT_ERROR
        if (size(y) < size(beta)) error stop FS_UNDERDEFINED_PROBLEM_ERROR
        if (.not.all(ieee_is_finite(x)) .or. &
            .not.all(ieee_is_finite(y)) .or. &
            .not.all(ieee_is_finite(beta))) error stop FS_INVALID_INPUT_ERROR
        significance = 0.05d0
        if (present(alpha)) significance = alpha
        if (.not.ieee_is_finite(significance)) error stop FS_INVALID_INPUT_ERROR
        if (significance <= 0.0d0 .or. significance >= 1.0d0) &
            error stop FS_INVALID_INPUT_ERROR
        if (present(options)) opts = options
        if (present(cov)) then
            if (size(cov,1) /= size(beta) .or. size(cov,2) /= size(beta)) &
                error stop FS_MATRIX_SIZE_ERROR
            cov = ieee_value(0.0d0, ieee_quiet_nan)
        end if
        if (present(stats) .or. present(cov)) then
            if (size(y) <= size(beta)) then
                if (opts%use_robust_weighting .or. method == FS_GLM_GAMMA) &
                    error stop FS_UNDERDEFINED_PROBLEM_ERROR
            end if
        end if
        if (present(stats)) then
            if (size(stats) /= size(beta)) error stop FS_ARRAY_SIZE_ERROR
            nan = ieee_value(0.0d0, ieee_quiet_nan)
            stats%standard_error = nan
            stats%wald_statistic = nan
            stats%p_value = nan
            stats%confidence_interval_lower = nan
            stats%confidence_interval_upper = nan
        end if

        ! Initialization
        select case (method)
        case (FS_GLM_BINOMIAL)
            if (any(y < 0.0d0) .or. any(y > 1.0d0)) &
                error stop FS_INVALID_INPUT_ERROR
            link_inv => glm_binomial_link_inv
        case (FS_GLM_POISSON)
            if (any(y < 0.0d0)) error stop FS_INVALID_INPUT_ERROR
            link_inv => glm_poisson_link_inv
        case (FS_GLM_GAMMA)
            if (any(y <= 0.0d0)) error stop FS_INVALID_INPUT_ERROR
            link_inv => glm_gamma_link_inv
        case default
            error stop FS_INVALID_INPUT_ERROR
        end select

        ! Process
        call irls_driver(method, x, y, beta, link_inv, weights = weights, &
            options = opts, info = fit_info, status = status, resid = resid)
        if (present(info)) info = fit_info
        if (present(stats) .or. present(cov)) then
            if (fit_info%converged) call coefficient_statistics(method, x, y, beta, &
            link_inv, opts, weights, significance, stats, cov)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    subroutine coefficient_statistics(method, x, y, beta, link_inv, opts, &
        weights, alpha, stats, cov)
        !! Recompute final-model covariance using scaled QR. Ordinary covariance
        !! is phi * (X**T W X)**(-1). Robust covariance is N/(N-P) times
        !! A**(-1) B A**(-T), where A differentiates the actual weighted score
        !! and B sums observation score outer products. Common score scaling
        !! cancels from this sandwich; column scaling is restored in log space.
        integer(int32), intent(in) :: method
            !! Validated family identifier.
        real(real64), intent(in), dimension(:,:) :: x
            !! Finite full-rank N-by-P design.
        real(real64), intent(in), dimension(:) :: y, beta
            !! Responses and converged coefficients.
        procedure(inverse_link_function), pointer, intent(in) :: link_inv
            !! Selected inverse-link callback.
        type(glm_options), intent(in) :: opts
            !! Validated weighting controls.
        procedure(robust_weight_function), pointer, intent(in), optional :: weights
            !! Optional deterministic observation-wise robust callback.
        real(real64), intent(in) :: alpha
            !! Validated two-sided significance level.
        type(glm_coefficient_statistics), intent(out), optional, dimension(:) :: stats
            !! P coefficient inference objects, on the linear-predictor scale.
        real(real64), intent(out), optional, dimension(:,:) :: cov
            !! P-by-P covariance, restoring normalization, column scales,
            !! and QR pivoting. Calculated without explicitly forming X**T W X.
        real(real64), allocatable, dimension(:) :: eta, mu, scale, column_scale, &
            logw, tau, se, wr, score, plus_score, minus_score, trial_eta, trial_mu
        real(real64), allocatable, dimension(:,:) :: design, matrix, qr, rhs, influence
        integer(int32), allocatable, dimension(:) :: pivots
        integer(int32) :: n, p, col, row, attempt, other
        real(real64) :: norm, logdisp, lognorm, h, cutoff, magnitude, critical, halfwidth
        real(real64) :: correlation, entry, left_norm, right_norm
        logical :: valid, use_t
        type(normal_distribution) :: normal
        type(t_distribution) :: student

        n = size(y)
        p = size(beta)
        allocate(eta(n), scale(p), column_scale(p), logw(n), se(p), pivots(p))
        call compute_predictor(x, beta, eta, valid)
        if (.not.valid) error stop FS_INVALID_INPUT_ERROR
        mu = link_inv(eta)
        scale = maxval(abs(x), dim = 1)
        design = x / spread(scale, dim = 1, ncopies = n)
        logdisp = 0.0d0
        use_t = method == FS_GLM_GAMMA .and. .not.opts%use_robust_weighting
        if (opts%use_robust_weighting) then
            call inference_weights(y - mu, opts, weights, wr)
            if (count(wr > 0.0d0) <= p) error stop FS_UNDERDEFINED_PROBLEM_ERROR
            lognorm = log(max(maxval(y), maxval(mu))) + log(maxval(wr))
            score = scaled_score(y - mu, wr, lognorm)
            allocate(matrix(p,p), trial_eta(n))
            do col = 1, p
                h = epsilon(1.0d0)**(1.0d0 / 3.0d0) * max(1.0d0, maxval(abs(eta)))
                do attempt = 1, 32
                    do row = 1, n
                        trial_eta(row) = bounded_sum(eta(row), h * design(row,col))
                    end do
                    trial_mu = link_inv(trial_eta)
                    valid = valid_mean(method, trial_mu)
                    if (valid) then
                        call inference_weights(y - trial_mu, opts, weights, wr)
                        plus_score = scaled_score(y - trial_mu, wr, lognorm)
                        do row = 1, n
                            trial_eta(row) = bounded_sum(eta(row), -h * design(row,col))
                        end do
                        trial_mu = link_inv(trial_eta)
                        valid = valid_mean(method, trial_mu)
                        if (valid) then
                            call inference_weights(y - trial_mu, opts, weights, wr)
                            minus_score = scaled_score(y - trial_mu, wr, lognorm)
                            exit
                        end if
                    end if
                    h = 0.5d0 * h
                end do
                if (.not.valid) error stop FS_INVALID_INPUT_ERROR
                matrix(:,col) = matmul(transpose(design), &
                    (plus_score - minus_score) / (2.0d0 * h))
            end do
            rhs = transpose(design * spread(score, dim = 2, ncopies = p))
            logdisp = log(real(n, real64) / real(n - p, real64))
        else
            select case (method)
            case (FS_GLM_BINOMIAL)
                logw = 0.5d0 * (log(mu) + log(1.0d0 - mu))
            case (FS_GLM_POISSON)
                logw = 0.5d0 * log(mu)
            case (FS_GLM_GAMMA)
                logw = log(mu)
                score = abs(y - mu)
                lognorm = -huge(1.0d0)
                do row = 1, n
                    if (score(row) > 0.0d0) lognorm = &
                        max(lognorm, log(score(row)) - log(mu(row)))
                end do
                if (lognorm == -huge(1.0d0)) then
                    logdisp = -huge(1.0d0)
                else
                    do row = 1, n
                        if (score(row) > 0.0d0) score(row) = &
                            bounded_exp(log(score(row)) - log(mu(row)) - lognorm)
                    end do
                    logdisp = 2.0d0 * (lognorm + log(norm2(score))) - log(real(n - p, real64))
                end if
            end select
            norm = maxval(logw)
            matrix = design * spread(bounded_exp(logw - norm), dim = 2, ncopies = p)
            allocate(rhs(p,p), source = 0.0d0)
            do col = 1, p
                rhs(col,col) = 1.0d0
            end do
            logdisp = logdisp - 2.0d0 * norm
        end if
        if (.not.all(ieee_is_finite(matrix))) error stop FS_INVALID_INPUT_ERROR
        column_scale = maxval(abs(matrix), dim = 1)
        if (any(column_scale == 0.0d0)) error stop FS_UNDERDEFINED_PROBLEM_ERROR
        matrix = matrix / spread(column_scale, dim = 1, ncopies = size(matrix, 1))
        pivots = 0
        call qr_factor(matrix, pivots, tau = tau, qr = qr)
        if (.not.all(ieee_is_finite(qr))) error stop FS_INVALID_INPUT_ERROR
        cutoff = epsilon(1.0d0) * real(max(n,p), real64) * abs(qr(1,1))
        do col = 1, p
            if (abs(qr(col,col)) <= cutoff) error stop FS_UNDERDEFINED_PROBLEM_ERROR
        end do
        if (opts%use_robust_weighting) rhs = mult_qr(.true., .true., qr, tau, rhs)
        allocate(influence(p,size(rhs,2)))
        do col = 1, size(rhs,2)
            influence(:,col) = solve_triangular_system(.true., .false., .true., qr(1:p,1:p), rhs(:,col))
        end do
        if (.not.all(ieee_is_finite(influence))) error stop FS_INVALID_INPUT_ERROR
        do col = 1, p
            row = pivots(col)
            magnitude = norm2(influence(col,:))
            se(row) = 0.0d0
            if (magnitude > 0.0d0) se(row) = bounded_exp(log(magnitude) + &
                0.5d0 * logdisp - log(scale(row)) - log(column_scale(row)))
            if (magnitude > 0.0d0 .and. logdisp > -0.5d0 * huge(1.0d0)) then
                if (se(row) == 0.0d0) error stop FS_INVALID_INPUT_ERROR
            end if
        end do
        if (.not.all(ieee_is_finite(se))) error stop FS_INVALID_INPUT_ERROR
        if (present(cov)) then
            cov = 0.0d0
            do col = 1, p
                left_norm = norm2(influence(col,:))
                if (left_norm == 0.0d0 .or. se(pivots(col)) == 0.0d0) cycle
                do other = col, p
                    right_norm = norm2(influence(other,:))
                    if (right_norm == 0.0d0 .or. se(pivots(other)) == 0.0d0) cycle
                    correlation = dot_product(influence(col,:) / left_norm, influence(other,:) / right_norm)
                    if (col == other) correlation = 1.0d0
                    if (correlation == 0.0d0) cycle
                    entry = sign(bounded_exp(log(abs(correlation)) + &
                        log(se(pivots(col))) + log(se(pivots(other)))), correlation)
                    if (.not.ieee_is_finite(entry)) error stop FS_INVALID_INPUT_ERROR
                    cov(pivots(col),pivots(other)) = entry
                    cov(pivots(other),pivots(col)) = entry
                end do
            end do
        end if
        if (.not.present(stats)) return
        normal%mean_value = 0.0d0
        normal%standard_deviation = 1.0d0
        student%dof = real(n - p, real64)
        critical = confidence_interval(normal, alpha, 1.0d0, 1)
        if (use_t) critical = confidence_interval(student, alpha, 1.0d0, 1)
        if (.not.ieee_is_finite(critical)) error stop FS_INVALID_INPUT_ERROR
        do col = 1, p
            stats(col)%standard_error = se(col)
            if (se(col) == 0.0d0) then
                if (beta(col) == 0.0d0) then
                    stats(col)%wald_statistic = 0.0d0
                else
                    stats(col)%wald_statistic = sign(ieee_value(0.0d0, ieee_positive_inf), beta(col))
                end if
            else if (beta(col) == 0.0d0) then
                stats(col)%wald_statistic = 0.0d0
            else
                magnitude = bounded_exp(log(abs(beta(col))) - log(se(col)))
                if (.not.ieee_is_finite(magnitude)) magnitude = ieee_value(0.0d0, ieee_positive_inf)
                stats(col)%wald_statistic = sign(magnitude, beta(col))
            end if
            stats(col)%p_value = 2.0d0 * normal%survival(abs(stats(col)%wald_statistic))
            if (use_t) stats(col)%p_value = 2.0d0 * student%survival(abs(stats(col)%wald_statistic))
            halfwidth = 0.0d0
            if (se(col) > 0.0d0) halfwidth = bounded_exp(log(critical) + log(se(col)))
            stats(col)%confidence_interval_lower = bounded_sum(beta(col), -halfwidth)
            stats(col)%confidence_interval_upper = bounded_sum(beta(col), halfwidth)
            if (.not.ieee_is_finite(stats(col)%confidence_interval_lower) .or. &
                .not.ieee_is_finite(stats(col)%confidence_interval_upper)) error stop FS_INVALID_INPUT_ERROR
        end do
    end subroutine

! ------------------------------------------------------------------------------
    subroutine inference_weights(r, opts, weights, wr)
        !! Evaluate and validate robust weights for inference perturbations.
        real(real64), intent(in), dimension(:) :: r
            !! Raw residuals at the proposed inference evaluation point.
        type(glm_options), intent(in) :: opts
            !! Validated robust controls.
        procedure(robust_weight_function), pointer, intent(in), optional :: weights
            !! Optional robust callback, defaulting to Tukey's biweight.
        real(real64), allocatable, intent(out), dimension(:) :: wr
            !! Allocated finite nonnegative weights, one per response.
        if (present(weights)) then
            if (associated(weights)) wr = weights(r, opts%robust_weighting_constant)
        end if
        if (.not.allocated(wr)) then
            if (present(weights)) then
                if (associated(weights)) error stop FS_INVALID_INPUT_ERROR
            end if
            wr = glm_tukey_weight(r, opts%robust_weighting_constant)
        end if
        if (size(wr) /= size(r)) error stop FS_ARRAY_SIZE_ERROR
        if (.not.all(ieee_is_finite(wr))) error stop FS_INVALID_INPUT_ERROR
        if (any(wr < 0.0d0)) error stop FS_INVALID_INPUT_ERROR
    end subroutine

! ------------------------------------------------------------------------------
    function scaled_score(r, wr, lognorm) result(score)
        !! Scale canonical-link observation scores by a common logarithmic
        !! factor. The Gamma sign cancels between bread and meat and is omitted.
        real(real64), intent(in), dimension(:) :: r, wr
            !! Raw residuals and finite nonnegative robust weights.
        real(real64), intent(in) :: lognorm
            !! Common log normalization, fixed throughout differentiation.
        real(real64), allocatable, dimension(:) :: score
            !! Scaled signed observation contributions, excluding the design.
        integer(int32) :: row
        allocate(score(size(r)), source = 0.0d0)
        do row = 1, size(r)
            if (wr(row) == 0.0d0 .or. r(row) == 0.0d0) cycle
            score(row) = sign(bounded_exp(log(abs(r(row))) + log(wr(row)) - lognorm), r(row))
        end do
        if (.not.all(ieee_is_finite(score))) error stop FS_INVALID_INPUT_ERROR
    end function

! ------------------------------------------------------------------------------
    subroutine irls_driver(method, x, y, beta, link_inv, weights, options, &
        info, status, resid)
        !! Solve the scaled weighted increment system. Family-specific
        !! logarithms avoid explicitly squaring derivatives or variances.
        integer(int32), intent(in) :: method
            !! Validated family identifier.
        real(real64), intent(in), dimension(:,:) :: x
            !! Validated N-by-P design matrix.
        real(real64), intent(in), dimension(:) :: y
            !! Validated N-element responses.
        real(real64), intent(inout), dimension(:) :: beta
            !! P initial coefficients, updated after each valid step.
        procedure(inverse_link_function), pointer, intent(in) :: link_inv
            !! Associated inverse-link callback for the selected family.
        procedure(robust_weight_function), pointer, intent(in), optional :: weights
            !! Optional robust callback; absent/null selects Tukey weighting.
        type(glm_options), intent(in), optional :: options
            !! Iteration controls, validated before allocation.
        type(glm_convergence_info), intent(out), optional :: info
            !! Completed-step diagnostics, reset before iteration.
        procedure(glm_iteration_update), intent(in), pointer, optional :: status
            !! An optional pointer to a routine that can be used to extract
            !! iteration information.
        real(real64), intent(out), optional, dimension(:) :: resid
            !! An optional P-element array where the residual will be written.

        ! Local Variables
        integer(int32) :: n, p, k, row, col, halvings
        integer(int32), allocatable, dimension(:) :: pivots
        real(real64), allocatable, dimension(:) :: eta, mu, trial_mu, logw, loggp, wr, &
            r, rhs, delta, beta_new, scale, wscale, tau, qrhs, solution
        real(real64), allocatable, dimension(:,:) :: wx, qr
        real(real64) :: diff, nan, wnorm, bnorm, magnitude, rank_tol
        logical :: valid
        type(glm_options) :: opts

        ! Initialization
        n = size(y)
        p = size(beta)
        nan = ieee_value(0.0d0, ieee_quiet_nan)
        if (present(options)) then
            opts = options
        else
            call opts%set_to_default()
        end if

        if (opts%max_iteration_count < 1) error stop FS_INVALID_INPUT_ERROR
        if (.not.ieee_is_finite(opts%tolerance)) &
            error stop FS_INVALID_INPUT_ERROR
        if (opts%tolerance <= 0.0d0) error stop FS_INVALID_INPUT_ERROR
        if (opts%use_robust_weighting) then
            if (.not.ieee_is_finite(opts%robust_weighting_constant)) &
                error stop FS_INVALID_INPUT_ERROR
            if (opts%robust_weighting_constant <= 0.0d0) &
                error stop FS_INVALID_INPUT_ERROR
        end if
        if (present(info)) then
            info%converged = .false.
            info%iteration_count = 0
            info%residual_value = nan
        end if
        if (present(resid)) then
            if (size(resid) /= p) error stop FS_ARRAY_SIZE_ERROR
        end if

        allocate(eta(n), logw(n), loggp(n), r(n), rhs(n), delta(p), &
            beta_new(p), scale(p), wscale(p), wx(n,p), pivots(p))
        scale = maxval(abs(x), dim = 1)
        if (any(scale == 0.0d0)) error stop FS_UNDERDEFINED_PROBLEM_ERROR

        ! Process
        do k = 1, opts%max_iteration_count
            call compute_predictor(x, beta, eta, valid)
            if (.not.valid) error stop FS_INVALID_INPUT_ERROR
            mu = link_inv(eta)
            if (.not.valid_mean(method, mu)) error stop FS_INVALID_INPUT_ERROR
            select case (method)
            case (FS_GLM_BINOMIAL)
                loggp = -log(mu) - log(1.0d0 - mu)
                logw = -0.5d0 * loggp
            case (FS_GLM_POISSON)
                loggp = -log(mu)
                logw = 0.5d0 * log(mu)
            case (FS_GLM_GAMMA)
                loggp = -2.0d0 * log(mu)
                logw = log(mu)
            end select
            r = y - mu
            if (opts%use_robust_weighting) then
                if (present(weights)) then
                    if (associated(weights)) then
                        wr = weights(r, opts%robust_weighting_constant)
                    else
                        wr = glm_tukey_weight(r, &
                            opts%robust_weighting_constant)
                    end if
                else
                    wr = glm_tukey_weight(r, &
                        opts%robust_weighting_constant)
                end if
                if (.not.allocated(wr)) error stop FS_INVALID_INPUT_ERROR
                if (size(wr) /= n) error stop FS_INVALID_INPUT_ERROR
                if (.not.all(ieee_is_finite(wr))) &
                    error stop FS_INVALID_INPUT_ERROR
                if (any(wr < 0.0d0)) error stop FS_INVALID_INPUT_ERROR
            else
                wr = spread(1.0d0, dim = 1, ncopies = n)
            end if

            if (count(wr > 0.0d0) < p) &
                error stop FS_UNDERDEFINED_PROBLEM_ERROR
            do row = 1, n
                if (wr(row) > 0.0d0) then
                    logw(row) = logw(row) + 0.5d0 * log(wr(row))
                else
                    logw(row) = -huge(1.0d0)
                end if
            end do

            wnorm = maxval(logw)
            rhs = -huge(1.0d0)
            do row = 1, n
                if (wr(row) == 0.0d0) cycle
                if (r(row) /= 0.0d0) rhs(row) = &
                    log(abs(r(row))) + loggp(row) + logw(row) - wnorm
            end do

            bnorm = max(0.0d0, maxval(rhs))
            do row = 1, n
                magnitude = bounded_exp(logw(row) - wnorm)
                wx(row,:) = (x(row,:) / scale) * magnitude
                rhs(row) = sign(bounded_exp(rhs(row) - bnorm), r(row))
            end do

            if (method == FS_GLM_GAMMA) rhs = -rhs
            wscale = maxval(abs(wx), dim = 1)
            if (any(wscale == 0.0d0)) error stop FS_UNDERDEFINED_PROBLEM_ERROR
            do col = 1, p
                wx(:,col) = wx(:,col) / wscale(col)
            end do
            pivots = 0
            call qr_factor(wx, pivots, tau = tau, qr = qr)
            if (.not.all(ieee_is_finite(qr))) error stop FS_INVALID_INPUT_ERROR
            rank_tol = epsilon(1.0d0) * real(max(n, p), real64) * abs(qr(1,1))
            do col = 1, p
                if (abs(qr(col,col)) <= rank_tol) &
                    error stop FS_UNDERDEFINED_PROBLEM_ERROR
            end do
            qrhs = mult_qr(.true., qr, tau, rhs)
            solution = solve_triangular_system(.true., .false., .true., &
                qr(1:p,1:p), qrhs(1:p))
            if (.not.all(ieee_is_finite(solution))) &
                error stop FS_INVALID_INPUT_ERROR

            delta(pivots) = solution
            do col = 1, p
                if (delta(col) == 0.0d0) cycle
                magnitude = log(abs(delta(col))) + bnorm - &
                    log(scale(col)) - log(wscale(col))
                delta(col) = sign(bounded_exp(magnitude), delta(col))
            end do

            if (.not.all(ieee_is_finite(delta))) error stop FS_INVALID_INPUT_ERROR
            do halvings = 0, 64
                do col = 1, p
                    beta_new(col) = bounded_sum(beta(col), delta(col))
                end do
                call compute_predictor(x, beta_new, eta, valid)
                if (valid) then
                    trial_mu = link_inv(eta)
                    valid = valid_mean(method, trial_mu)
                    if (valid) valid = acceptable_step(method, y, mu, trial_mu, wr)
                end if
                if (valid) exit
                delta = 0.5d0 * delta
            end do
            if (.not.valid) error stop FS_INVALID_INPUT_ERROR

            ! Convergence check
            diff = maxval(abs(delta))
            beta = beta_new
            if (present(info)) then
                info%residual_value = diff
                info%iteration_count = k
            end if
            if (diff < opts%tolerance .and. halvings == 0) then
                if (present(info)) then
                    info%converged = .true.
                end if
                exit
            end if
            if (present(status)) then
                if (associated(status)) then
                    call status(k, delta, beta)
                end if
            end if
        end do

        ! Return additional outputs
        if (present(resid)) resid = delta
    end subroutine

! ------------------------------------------------------------------------------
    pure elemental function bounded_exp(value) result(rst)
        !! Exponentiate without overflow; retain representable subnormals and
        !! flush smaller values to zero. Return
        !! NaN for a nonfinite argument or an unrepresentably large result.
        real(real64), intent(in) :: value
            !! Logarithm of a proposed nonnegative value.
        real(real64) :: rst
            !! Representable exponential, zero on underflow, or quiet NaN.
        rst = ieee_value(0.0d0, ieee_quiet_nan)
        if (.not.ieee_is_finite(value)) return
        if (value >= log(huge(1.0d0))) return
        rst = 0.0d0
        if (value >= log(tiny(1.0d0)) + log(epsilon(1.0d0))) rst = exp(value)
    end function

! ------------------------------------------------------------------------------
    pure function acceptable_step(method, y, mu, trial_mu, wr) result(valid)
        !! Compare half-deviances on a common scale, excluding common
        !! dispersion. Robust weights are fixed during backtracking. A small
        !! roundoff allowance prevents rejecting steps at a numerical optimum.
        integer(int32), intent(in) :: method
            !! Validated family identifier.
        real(real64), intent(in), dimension(:) :: y, mu, trial_mu, wr
            !! Responses, current/trial valid means, and finite robust weights.
        logical :: valid
            !! True if the proposed fixed-weight deviance does not increase
            !! beyond 32 * epsilon * max(1, current scaled half-deviance).
        real(real64) :: old_cost, new_cost, norm, lognorm, weight, old_term, new_term
        real(real64) :: logratio, logtrial, factor, weight_norm
        integer(int32) :: row
        old_cost = 0.0d0
        new_cost = 0.0d0
        norm = 1.0d0
        lognorm = 0.0d0
        weight_norm = maxval(wr)
        if (method == FS_GLM_POISSON) then
            norm = max(maxval(y, mask = wr > 0.0d0), &
                maxval(mu, mask = wr > 0.0d0), maxval(trial_mu, mask = wr > 0.0d0))
        else if (method == FS_GLM_GAMMA) then
            lognorm = max(0.0d0, maxval(log(y) - log(mu), mask = wr > 0.0d0), &
                maxval(log(y) - log(trial_mu), mask = wr > 0.0d0))
        end if
        factor = bounded_exp(-lognorm)
        do row = 1, size(y)
            if (wr(row) == 0.0d0) cycle
            weight = wr(row) / weight_norm
            select case (method)
            case (FS_GLM_BINOMIAL)
                old_term = -y(row) * log(mu(row)) - &
                    (1.0d0 - y(row)) * log(1.0d0 - mu(row))
                new_term = -y(row) * log(trial_mu(row)) - &
                    (1.0d0 - y(row)) * log(1.0d0 - trial_mu(row))
            case (FS_GLM_POISSON)
                old_term = mu(row) / norm
                new_term = trial_mu(row) / norm
                if (y(row) > 0.0d0) then
                    old_term = old_term + (y(row) / norm) * &
                        (log(y(row)) - log(mu(row)) - 1.0d0)
                    new_term = new_term + (y(row) / norm) * &
                        (log(y(row)) - log(trial_mu(row)) - 1.0d0)
                end if
            case (FS_GLM_GAMMA)
                logratio = log(y(row)) - log(mu(row))
                logtrial = log(y(row)) - log(trial_mu(row))
                old_term = bounded_exp(logratio - lognorm) - factor * (1.0d0 + logratio)
                new_term = bounded_exp(logtrial - lognorm) - factor * (1.0d0 + logtrial)
            end select
            old_cost = old_cost + weight * old_term
            new_cost = new_cost + weight * new_term
        end do
        valid = new_cost <= old_cost + &
            32.0d0 * epsilon(1.0d0) * max(1.0d0, abs(old_cost))
    end function

! ------------------------------------------------------------------------------
    pure function bounded_sum(left, right) result(rst)
        !! Add finite operands only when their sum is representable.
        real(real64), intent(in) :: left, right
            !! Proposed summands.
        real(real64) :: rst
            !! Sum, or quiet NaN if nonfinite inputs/overflow are detected.
        rst = ieee_value(0.0d0, ieee_quiet_nan)
        if (.not.ieee_is_finite(left) .or. .not.ieee_is_finite(right)) return
        if (right > 0.0d0) then
            if (left > huge(1.0d0) - right) return
        else if (right < 0.0d0) then
            if (left < -huge(1.0d0) - right) return
        end if
        rst = left + right
    end function

! ------------------------------------------------------------------------------
    pure subroutine compute_predictor(x, beta, eta, valid)
        !! Compute X beta with checked products and sums before overflow.
        !! Cancellation requiring unrepresentable intermediate terms is rejected.
        real(real64), intent(in), dimension(:,:) :: x
            !! Finite N-by-P design.
        real(real64), intent(in), dimension(:) :: beta
            !! P proposed coefficients.
        real(real64), intent(out), dimension(:) :: eta
            !! N predictors; only usable if valid is true.
        logical, intent(out) :: valid
            !! True when all arithmetic produces finite predictors.
        integer(int32) :: row, col
        real(real64) :: term
        valid = .false.
        if (.not.all(ieee_is_finite(beta))) return
        eta = 0.0d0
        do col = 1, size(beta)
            do row = 1, size(eta)
                if (abs(x(row,col)) > 1.0d0) then
                    if (abs(beta(col)) > huge(1.0d0) / abs(x(row,col))) return
                end if
                term = x(row,col) * beta(col)
                eta(row) = bounded_sum(eta(row), term)
                if (.not.ieee_is_finite(eta(row))) return
            end do
        end do
        valid = .true.
    end subroutine

! ------------------------------------------------------------------------------
    pure function valid_mean(method, mu) result(valid)
        !! Reject unrepresentable means and rounded binomial boundary values.
        integer(int32), intent(in) :: method
            !! Validated family identifier.
        real(real64), intent(in), dimension(:) :: mu
            !! N inverse-link results.
        logical :: valid
            !! True for finite positive means, also < 1 for binomial.
        valid = all(ieee_is_finite(mu))
        if (.not.valid) return
        valid = all(mu > 0.0d0)
        if (method == FS_GLM_BINOMIAL) valid = valid .and. all(mu < 1.0d0)
    end function

! ******************************************************************************
! FS_GLM_BINOMIAL ROUTINES
! ------------------------------------------------------------------------------
    pure function glm_binomial_link_inv(eta) result(mu)
        !! Evaluate the logistic inverse link without exponent overflow.
        real(real64), intent(in), dimension(:) :: eta
            !! N finite log-odds.
        real(real64), allocatable, dimension(:) :: mu
            !! N probabilities; the driver rejects rounded zero/one values.
        integer(int32) :: row
        real(real64) :: term
        allocate(mu(size(eta)))
        do row = 1, size(eta)
            term = bounded_exp(-abs(eta(row)))
            if (eta(row) >= 0.0d0) then
                mu(row) = 1.0d0 / (1.0d0 + term)
            else
                mu(row) = term / (1.0d0 + term)
            end if
        end do
    end function

! ******************************************************************************
! FS_GLM_POISSON ROUTINES
! ------------------------------------------------------------------------------
    pure function glm_poisson_link_inv(eta) result(mu)
        !! Invert the Poisson log link with bounded exponentiation.
        real(real64), intent(in), dimension(:) :: eta
            !! N finite log-means.
        real(real64), allocatable, dimension(:) :: mu
            !! N means exp(eta); overflow is NaN and underflow is zero.
        mu = bounded_exp(eta)
    end function

! ******************************************************************************
! FS_GLM_GAMMA ROUTINES
! ------------------------------------------------------------------------------
    pure function glm_gamma_link_inv(eta) result(mu)
        !! Invert the reciprocal Gamma link g(mu) = 1 / mu.
        real(real64), intent(in), dimension(:) :: eta
            !! N predictors; nonpositive or unrepresentable values yield NaN.
        real(real64), allocatable, dimension(:) :: mu
            !! N means 1 / eta.
        integer(int32) :: row
        allocate(mu(size(eta)))
        mu = ieee_value(0.0d0, ieee_quiet_nan)
        do row = 1, size(eta)
            if (eta(row) > 1.0d0 / huge(1.0d0)) mu(row) = 1.0d0 / eta(row)
        end do
    end function

! ******************************************************************************
! DEFAULT ROBUST WEIGHTING ROUTINES
! ------------------------------------------------------------------------------
    pure function glm_tukey_weight(u, k) result(rst)
        !! Adapt Tukey's elemental biweight to the vector callback interface:
        !! (1 - (u/k)**2)**2 for abs(u) <= k, and zero otherwise.
        real(real64), intent(in), dimension(:) :: u
            !! N raw response residuals, without scale standardization.
        real(real64), intent(in) :: k
            !! Positive cutoff in response units.
        real(real64), allocatable, dimension(:) :: rst
            !! N nonnegative multiplicative weights.
        rst = tukey_biweight_irls_weight(u, k)
    end function

! ------------------------------------------------------------------------------
end module
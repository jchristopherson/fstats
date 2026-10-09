module fstats_glm_irls
    !! Numerically scaled IRLS for logit-binomial, log-Poisson, and
    !! reciprocal-Gamma generalized linear models.
    use iso_fortran_env
    use ieee_arithmetic
    use fstats_glm_types
    use fstats_errors
    use fstats_robust_statistics
    use linalg, only : qr_factor, mult_qr, solve_triangular_system
    implicit none
    private
    public :: irls

contains
! ------------------------------------------------------------------------------
    subroutine irls(method, x, y, beta, weights, options, info)
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

        ! Local Variables
        procedure(inverse_link_function), pointer :: link_inv

        ! Input Validation
        if (size(x, 1) /= size(y) .or. size(x, 2) /= size(beta)) &
            error stop FS_MATRIX_SIZE_ERROR
        if (size(y) == 0 .or. size(beta) == 0) &
            error stop FS_INVALID_INPUT_ERROR
        if (size(y) < size(beta)) error stop FS_UNDERDEFINED_PROBLEM_ERROR
        if (.not.all(ieee_is_finite(x)) .or. &
            .not.all(ieee_is_finite(y)) .or. &
            .not.all(ieee_is_finite(beta))) error stop FS_INVALID_INPUT_ERROR

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
        call irls_driver(method, x, y, beta, link_inv, &
            weights = weights, options = options, info = info)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine irls_driver(method, x, y, beta, link_inv, weights, options, info)
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
        end do
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
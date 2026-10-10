program example
    use iso_fortran_env
    use fstats
    use fplot_core
    implicit none

    ! Local Variables
    integer(int32), parameter :: eval_size = 500
    integer(int32) :: i
    real(real64) :: dx, x(10,2), y(10), beta_p(2), eta(eval_size), &
        f(eval_size), xi(eval_size,2)
    type(glm_coefficient_statistics) :: stats(2)
    type(glm_convergence_info) :: info
    type(glm_options) :: opts

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pdata

    ! The model is log(p / (1 - p)) = beta_p(1) + beta_p(2) * x,
    ! where p = P(y = 1). The first design column supplies the intercept.
    dx = 0.5d0
    do i = 1, size(x, 1)
        x(i,1) = 1.0d0
        x(i,2) = (i - 1.0d0) * dx
    end do
    y = [0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 1.0d0, 1.0d0, 0.0d0, 1.0d0, 1.0d0]
    beta_p = 0.0d0
    xi(:,1) = 1.0d0
    xi(:,2) = linspace(minval(x(:,2)), maxval(x(:,2)), eval_size)

    ! Disable robust weighting for an ordinary Bernoulli likelihood fit.
    opts%use_robust_weighting = .false.

    ! Solve the problem
    call irls(FS_GLM_BINOMIAL, x, y, beta_p, options = opts, stats = stats, &
        info = info)

    ! Converged means the maximum absolute coefficient update met the tolerance.
    ! Iteration count includes completed updates, not backtracking step halvings.
    ! residual_value is the last maximum coefficient change, not a y residual.
    print "(A, L0)", new_line('a') // "CONVERGENCE STATUS: ", info%converged
    print "(A, I0)", "ITERATION COUNT: ", info%iteration_count
    print "(A, ES12.5)", "MAX COEFFICIENT CHANGE: ", info%residual_value

    ! Coefficient 1 is the log-odds at x = 0; coefficient 2 is the log-odds
    ! change for a one-unit increase in x. exp(beta_p(2)) is the odds ratio.
    ! These values are on the linear-predictor scale, not probabilities.
    ! Statistics are NaN if the fit did not converge.
    print "(A)", new_line('a') // "MODEL OUTPUT:"
    do i = 1, size(beta_p)
        print "(A, I0)", "COEFFICIENT: ", i
        print "(A, F0.5)", achar(9) // "Value: ", beta_p(i)
        ! Standard error measures coefficient uncertainty from the fitted
        ! model's information covariance, assuming independent observations.
        print "(A, F0.5)", achar(9) // "Std. Error: ", stats(i)%standard_error
        ! Wald z = coefficient / standard error tests a zero coefficient
        ! using an asymptotic standard normal reference distribution.
        print "(A, F0.5)", achar(9) // "Wald Stat: ", stats(i)%wald_statistic
        ! These are approximate 95% Wald bounds (default alpha = 0.05),
        ! not bounds on the predicted probability or an interval half-width.
        print "(A, F0.5)", achar(9) // "Upper CI: ", stats(i)%confidence_interval_upper
        print "(A, F0.5)", achar(9) // "Lower CI: ", stats(i)%confidence_interval_lower
        ! The two-sided p-value is a null-hypothesis tail probability,
        ! not the probability that the coefficient is zero or unimportant.
        print "(A, F0.5)", achar(9) // "P-Value: ", stats(i)%p_value
    end do

    ! eta contains fitted log-odds; the inverse logit converts them to
    ! probabilities on a dense grid. The curve has no uncertainty band.
    eta = matmul(xi, beta_p)
    f = glm_binomial_link_inv(eta)

    ! Plot
    call plt%initialize()
    call plt%set_title("Binomial GLM with Logit Link")
    call plt%set_x_axis_title("Predictor x")
    call plt%set_y_axis_title("Response and P(y = 1)")

    call pdata%define_data(x(:,2), y)
    call pdata%set_draw_line(.false.)
    call pdata%set_draw_markers(.true.)
    call pdata%set_name("Observed binary response")
    call plt%push(pdata)

    call plt%push(xi(:,2), f, name = "Fitted P(y = 1)")

    call plt%draw()
end program
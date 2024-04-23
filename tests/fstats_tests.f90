program tests
    use iso_fortran_env
    use fstats_distribution_tests
    use fstats_statistics_tests
    use fstats_regression_tests
    use fstats_experimental_design_tests
    use fstats_nonlinear_regression_tests
    use fstats_allan_tests
    implicit none

    ! Variables
    logical :: local, overall

    ! Initialization
    overall = .true.

    ! Distribution Tests
    local = t_distribution_test_1()
    if (.not.local) overall = .false.

    local = normal_distribution_test_1()
    if (.not.local) overall = .false.

    local = f_distribution_test_1()
    if (.not.local) overall = .false.

    local = chi_squared_distribution_test_1()
    if (.not.local) overall = .false.

    local = binomial_distribution_test_1()
    if (.not.local) overall = .false.

    ! Statistics Tests
    local = mean_test_1()
    if (.not.local) overall = .false.

    local = variance_test_1()
    if (.not.local) overall = .false.

    local = standard_deviation_test_1()
    if (.not.local) overall = .false.

    local = median_test_1()
    if (.not.local) overall = .false.

    local = r_squared_test_1()
    if (.not.local) overall = .false.

    local = t_test_test_1()
    if (.not.local) overall = .false.

    local = f_test_test_1()
    if (.not.local) overall = .false.

    local = anova_test_1()
    if (.not.local) overall = .false.

    local = anova_test_2()
    if (.not.local) overall = .false.

    local = anova_test_3()
    if (.not.local) overall = .false.

    local = confidence_interval_test_1()
    if (.not.local) overall = .false.

    local = beta_test_1()
    if (.not.local) overall = .false.

    local = beta_test_2()
    if (.not.local) overall = .false.

    local = incomplete_gamma_test_1()
    if (.not.local) overall = .false.

    local = design_matrix_test_1()
    if (.not.local) overall = .false.

    local = regression_test_1()
    if (.not.local) overall = .false.

    local = bootstrap_regression_test_1()
    if (.not.local) overall = .false.

    ! Experimental Design
    ! local = get_full_matrix_size_test_1()
    ! if (.not.local) overall = .false.

    ! local = full_factorial_test_1()
    ! if (.not.local) overall = .false.

    ! ! Nonlinear Regression
    ! local = test_prototype_function_call()
    ! if (.not.local) overall = .false.

    ! local = test_jacobian()
    ! if (.not.local) overall = .false.

    ! local = test_nl_least_squares()
    ! if (.not.local) overall = .false.

    ! local = test_bootstrap_nl_least_squares()
    ! if (.not.local) overall = .false.

    ! ! Allan Variance Tests
    ! local = test_allan_variance()
    ! if (.not.local) overall = .false.

    ! ! Trimmed Mean Tests
    ! local = trimmed_mean_test_1()
    ! if (.not.local) overall = .false.

    ! ! Covariance Tests
    ! local = test_covariance_1()
    ! if (.not.local) overall = .false.

    ! ! Correlation Tests
    ! local = test_correlation_1()
    ! if (.not.local) overall = .false.

    ! ! Additional Tests
    ! local = test_pooled_variance_1()
    ! if (.not.local) overall = .false.

    ! local = test_bartlett_1()
    ! if (.not.local) overall = .false.

    ! local = test_levene_1()
    ! if (.not.local) overall = .false.

    ! local = test_standardized_variable()
    ! if (.not.local) overall = .false.

    ! local = test_sample_size()
    ! if (.not.local) overall = .false.

    ! End
    if (.not.overall) then
        stop 1
    end if
end program
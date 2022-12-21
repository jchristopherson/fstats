program tests
    use iso_fortran_env
    use fstats_distribution_tests
    use fstats_statistics_tests
    use fstats_regression_tests
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

    local = coefficient_matrix_test_1()
    if (.not.local) overall = .false.

    local = regression_test_1()
    if (.not.local) overall = .false.

    ! End
    if (overall) then
        call exit(0)
    else
        call exit(1)
    end if
end program
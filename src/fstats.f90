module fstats
    !! FSTATS is a modern Fortran statistical library containing routines for 
    !! computing basic statistical properties, hypothesis testing, regression, 
    !! special functions, and experimental design.
    use iso_fortran_env
    use fstats_special_functions
    use fstats_descriptive_statistics
    use fstats_hypothesis
    use fstats_distributions
    use fstats_anova
    use fstats_helper_routines
    use fstats_regression
    use fstats_experimental_design
    use fstats_allan
    use fstats_bootstrap
    use fstats_sampling
    use fstats_smoothing
    use fstats_mcmc
    ! use fstats_mcmc_fitting
    implicit none
    private

    ! FSTATS_DISTRIBUTION.F90
    public :: distribution
    public :: distribution_function
    public :: distribution_property
    public :: distribution_recenter
    public :: t_distribution
    public :: normal_distribution
    public :: f_distribution
    public :: chi_squared_distribution
    public :: binomial_distribution
    public :: multivariate_distribution
    public :: multivariate_distribution_function
    public :: multivariate_normal_distribution
    public :: log_normal_distribution
    public :: poisson_distribution

    ! FSTATS_DESCRIPTIVE_STATISTICS.F90
    public :: mean
    public :: variance
    public :: standard_deviation
    public :: median
    public :: quantile
    public :: trimmed_mean
    public :: covariance
    public :: pooled_variance

    ! FSTATS_SPECIAL_FUNCTIONS
    public :: beta
    public :: regularized_beta
    public :: incomplete_beta
    public :: incomplete_gamma_lower
    public :: incomplete_gamma_upper
    public :: digamma

    ! FSTATS_SMOOTHING.F90
    public :: lowess

    ! FSTATS_SAMPLING.F90
    public :: box_muller_sample
    public :: rejection_sample
    public :: sample_normal_multivariate
    
    ! FSTATS_MCMC.F90
    public :: chain_builder
    public :: mcmc_sampler
    public :: mcmc_target
    public :: evaluate_model
    public :: mcmc_proposal
   
    ! FSTATS_REGRESSION.F90
    public :: iteration_controls
    public :: convergence_info
    public :: lm_solver_options
    public :: regression_function
    public :: iteration_update
    public :: regression_statistics
    public :: r_squared
    public :: adjusted_r_squared
    public :: correlation
    public :: design_matrix
    public :: covariance_matrix
    public :: linear_least_squares
    public :: calculate_regression_statistics
    public :: jacobian
    public :: nonlinear_least_squares
    public :: FS_LEVENBERG_MARQUARDT_UPDATE
    public :: FS_QUADRATIC_UPDATE
    public :: FS_NIELSEN_UPDATE

    ! FSTATS_HYPOTHESIS.F90
    public :: confidence_interval
    public :: t_test_equal_variance
    public :: t_test_unequal_variance
    public :: t_test_paired
    public :: f_test
    public :: bartletts_test
    public :: levenes_test
    public :: sample_size

    ! FSTATS_HELPER_ROUTINES.F90
    public :: difference
    public :: factorial

    ! FSTATS_EXPERIMENTAL_DESIGN.F90
    public :: get_full_factorial_matrix_size
    public :: full_factorial
    public :: doe_fit_model
    public :: doe_evaluate_model
    public :: doe_model

    ! FSTATS_BOOTSTRAP.F90
    public :: bootstrap_resampling_routine
    public :: bootstrap_statistic_routine
    public :: random_resample
    public :: scaled_random_resample
    public :: bootstrap_statistics
    public :: bootstrap

    ! FSTATS_ANOVA.F90
    public :: anova_factor
    public :: single_factor_anova_table
    public :: two_factor_anova_table
    public :: anova

    ! FSTATS_ALLAN.F90
    public :: allan_variance
end module
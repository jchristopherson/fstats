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
    public :: t_distribution
    public :: normal_distribution
    public :: f_distribution
    public :: chi_squared_distribution
    public :: binomial_distribution
    public :: multivariate_distribution
    public :: multivariate_distribution_function
    public :: multivariate_normal_distribution
    public :: log_normal_distribution

    public :: mean
    public :: variance
    public :: standard_deviation
    public :: median
    public :: covariance
    public :: r_squared
    public :: adjusted_r_squared
    public :: correlation
    public :: quantile
    public :: t_test_equal_variance
    public :: t_test_unequal_variance
    public :: t_test_paired
    public :: f_test
    public :: anova
    public :: anova_factor
    public :: single_factor_anova_table
    public :: two_factor_anova_table
    public :: confidence_interval
    public :: beta
    public :: regularized_beta
    public :: incomplete_beta
    public :: digamma
    public :: incomplete_gamma_upper
    public :: incomplete_gamma_lower
    public :: design_matrix
    public :: covariance_matrix
    public :: linear_least_squares
    public :: regression_statistics
    public :: get_full_factorial_matrix_size
    public :: full_factorial
    public :: iteration_controls
    public :: lm_solver_options
    public :: convergence_info
    public :: regression_function
    public :: iteration_update
    public :: jacobian
    public :: nonlinear_least_squares
    public :: allan_variance
    public :: trimmed_mean
    public :: difference
    public :: factorial
    public :: bootstrap_resampling_routine
    public :: bootstrap_statistic_routine
    public :: random_resample
    public :: scaled_random_resample
    public :: bootstrap_statistics
    public :: bootstrap
    
    public :: lowess
    public :: pooled_variance
    public :: bartletts_test
    public :: levenes_test
    public :: sample_size
    public :: FS_LEVENBERG_MARQUARDT_UPDATE
    public :: FS_QUADRATIC_UPDATE
    public :: FS_NIELSEN_UPDATE
    public :: doe_fit_model
    public :: doe_evaluate_model
    public :: doe_model

    ! FSTATS_SAMPLING.F90
    public :: box_muller_sample
    public :: rejection_sample
    public :: sample_normal_multivariate

    ! FSTATS_MCMC.F90
    public :: metropolis_hastings

    ! FSTATS_MCMC_FITTING.F90
    ! public :: mcmc_regression
   
end module
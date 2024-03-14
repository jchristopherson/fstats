module fstats
    !! FSTATS is a modern Fortran statistical library containing routines for 
    !! computing basic statistical properties, hypothesis testing, regression, 
    !! special functions, and experimental design.
    use iso_fortran_env
    use fstats_special_functions
    use fstats_descriptive_statistics
    use fstats_distributions
    use fstats_anova
    use fstats_helper_routines
    use fstats_regression
    use fstats_experimental_design
    use fstats_allan
    use fstats_bootstrap
    use fstats_sampling
    use ferror
    implicit none
    private
    public :: distribution
    public :: distribution_function
    public :: distribution_property
    public :: t_distribution
    public :: normal_distribution
    public :: f_distribution
    public :: chi_squared_distribution
    public :: binomial_distribution
    public :: mean
    public :: variance
    public :: standard_deviation
    public :: median
    public :: r_squared
    public :: adjusted_r_squared
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
    public :: coefficient_matrix
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
    public :: bootstrap_regression_statistics
    public :: bootstrap_linear_least_squares
    public :: box_muller_sample
    public :: rejection_sample
    public :: FS_LEVENBERG_MARQUARDT_UPDATE
    public :: FS_QUADRATIC_UPDATE
    public :: FS_NIELSEN_UPDATE
   
end module
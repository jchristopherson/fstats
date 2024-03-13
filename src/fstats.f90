module fstats
    !! FSTATS is a modern Fortran statistical library containing routines for 
    !! computing basic statistical properties, hypothesis testing, regression, 
    !! special functions, and experimental design.
    use iso_fortran_env
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
    public :: bootstrap_regression_statistics
    public :: bootstrap_linear_least_squares
    public :: box_muller_sample
    public :: rejection_sample
    public :: FS_LEVENBERG_MARQUARDT_UPDATE
    public :: FS_QUADRATIC_UPDATE
    public :: FS_NIELSEN_UPDATE
    public :: assignment(=)


! ******************************************************************************
! CONSTANTS
! ------------------------------------------------------------------------------
    integer(int32), parameter :: FS_LEVENBERG_MARQUARDT_UPDATE = 1
    integer(int32), parameter :: FS_QUADRATIC_UPDATE = 2
    integer(int32), parameter :: FS_NIELSEN_UPDATE = 3

! ******************************************************************************
! OPERATORS
! ------------------------------------------------------------------------------
    interface assignment (=)
        module procedure :: ic_equal
        module procedure :: ci_equal
        module procedure :: lso_equal
    end interface
    
! ******************************************************************************
! DISTRIBUTIONS
! ------------------------------------------------------------------------------
    type, abstract :: distribution
        !! Defines a probability distribution.
    contains
        procedure(distribution_function), deferred, pass :: pdf
            !! Computes the probability density function.
        procedure(distribution_function), deferred, pass :: cdf
            !! Computes the cumulative distribution function.
        procedure(distribution_property), deferred, pass :: mean
            !! Computes the mean of the distribution.
        procedure(distribution_property), deferred, pass :: median
            !! Computes the median of the distribution.
        procedure(distribution_property), deferred, pass :: mode
            !! Computes the mode of the distribution.
        procedure(distribution_property), deferred, pass :: variance
            !! Computes the variance of the distribution.
    end type

    interface
        pure elemental function distribution_function(this, x) result(rst)
            !! Defines the interface for a probability distribution function.
            use iso_fortran_env, only : real64
            import distribution
            class(distribution), intent(in) :: this
                !! The distribution object.
            real(real64), intent(in) :: x
                !! The value at which to evaluate the function.
            real(real64) :: rst
                !! The value of the function.
        end function

        pure function distribution_property(this) result(rst)
            !! Computes the value of a distribution property.
            use iso_fortran_env, only : real64
            import distribution
            class(distribution), intent(in) :: this
                !! The distribution object.
            real(real64) :: rst
                !! The property value.
        end function
    end interface

    type, extends(distribution) :: t_distribution
        !! Defines Student's T-Distribution.
        real(real64) :: dof
            !! The number of degrees of freedom.
    contains
        procedure, public :: pdf => td_pdf
        procedure, public :: cdf => td_cdf
        procedure, public :: mean => td_mean
        procedure, public :: median => td_median
        procedure, public :: mode => td_mode
        procedure, public :: variance => td_variance
    end type

    ! distributions_t.f90
    interface
        pure module elemental function td_pdf(this, x) result(rst)
            !! Computes the probability density function.
            !!
            !! The PDF for Student's T-Distribution is given as 
            !! $$ f(t) = \frac{ \Gamma \left( \frac{\nu + 1}{2} \right) }
            !! { \sqrt{\nu \pi} \Gamma \left( \frac{\nu}{2} \right) } 
            !! \left( 1 + \frac{t^2}{\nu} \right)^{-(\nu + 1) / 2} $$
            class(t_distribution), intent(in) :: this
                !! The t_distribution object.
            real(real64), intent(in) :: x
                !! The value at which to evaluate the function.
            real(real64) :: rst
                !! The value of the function.
        end function

        pure module elemental function td_cdf(this, x) result(rst)
            !! Computes the cumulative distribution function.
            !!
            !! The CDF for Student's T-Distribution is given as
            !! $$ F(t) = \int_{-\infty}^{t} f(u) \,du = 1 - \frac{1}{2} I_{x(t)}
            !! \left( \frac{\nu}{2}, \frac{1}{2} \right) $$
            !! where $$ x(t) = \frac{\nu}{\nu + t^2} $$
            class(t_distribution), intent(in) :: this
                !! The t_distribution object.
            real(real64), intent(in) :: x
                !! The value at which to evaluate the function.
            real(real64) :: rst
                !! The value of the function.
        end function

        pure module function td_mean(this) result(rst)
            !! Computes the mean of the distribution.
            class(t_distribution), intent(in) :: this
                !! The t_distribution object.
            real(real64) :: rst
                !! The mean.
        end function

        pure module function td_median(this) result(rst)
            !! Computes the median of the distribution.
            class(t_distribution), intent(in) :: this
                !! The t_distribution object.
            real(real64) :: rst
                !! The median.
        end function

        pure module function td_mode(this) result(rst)
            !! Computes the mode of the distribution.
            class(t_distribution), intent(in) :: this
                !! The t_distribution object.
            real(real64) :: rst
                !! The mode.
        end function

        pure module function td_variance(this) result(rst)
            !! Computes the variance of the distribution.
            class(t_distribution), intent(in) :: this
                !! The t_distribution object.
            real(real64) :: rst
                !! The variance.
        end function
    end interface

    type, extends(distribution) :: normal_distribution
        !! Defines a normal distribution.
        real(real64) :: standard_deviation
            !! The standard deviation of the distribution.
        real(real64) :: mean_value
            !! The mean value of the distribution.
    contains
        procedure, public :: pdf => nd_pdf
        procedure, public :: cdf => nd_cdf
        procedure, public :: mean => nd_mean
        procedure, public :: median => nd_median
        procedure, public :: mode => nd_mode
        procedure, public :: variance => nd_variance
        procedure, public :: standardize => nd_standardize
    end type

    ! distributions_normal.f90
    interface
        pure module elemental function nd_pdf(this, x) result(rst)
            !! Computes the probability density function.
            !!
            !! The PDF for a normal distribution is given as 
            !! $$ f(x) = \frac{1}{\sigma \sqrt{2 \pi}} \exp \left(-\frac{1}{2} 
            !! \left( \frac{x - \mu}{\sigma} \right)^2 \right) $$.
            class(normal_distribution), intent(in) :: this
                !! The normal_distribution object.
            real(real64), intent(in) :: x
                !! The value at which to evaluate the function.
            real(real64) :: rst
                !! The value of the function.
        end function

        pure module elemental function nd_cdf(this, x) result(rst)
            !! Computes the cumulative distribution function.
            !!
            !! The CDF for a normal distribution is given as 
            !! $$ F(x) = \frac{1}{2} \left( 1 + erf \left( \frac{x - \mu}
            !! {\sigma \sqrt{2}} \right) \right) $$
            class(normal_distribution), intent(in) :: this
                !! The normal_distribution object.
            real(real64), intent(in) :: x
                !! The value at which to evaluate the function.
            real(real64) :: rst
                !! The value of the function.
        end function

        pure module function nd_mean(this) result(rst)
            !! Computes the mean of the distribution.
            class(normal_distribution), intent(in) :: this
                !! The normal_distribution object.
            real(real64) :: rst
                !! The mean.
        end function

        pure module function nd_median(this) result(rst)
            !! Computes the median of the distribution.
            class(normal_distribution), intent(in) :: this
                !! The normal_distribution object.
            real(real64) :: rst
                !! The median.
        end function

        pure module function nd_mode(this) result(rst)
            !! Computes the mode of the distribution.
            class(normal_distribution), intent(in) :: this
                !! The normal_distribution object.
            real(real64) :: rst
                !! The mode.
        end function

        pure module function nd_variance(this) result(rst)
            !! Computes the variance of the distribution.
            class(normal_distribution), intent(in) :: this
                !! The normal_distribution object.
            real(real64) :: rst
                !! The variance.
        end function

        module subroutine nd_standardize(this)
            !! Standardizes the normal distribution to a mean of 0 and a 
            !! standard deviation of 1.
            class(normal_distribution), intent(inout) :: this
                !! The normal_distribution object.
        end subroutine
    end interface

    type, extends(distribution) :: f_distribution
        !! Defines an F-distribution.
        real(real64) :: d1
            !! The measure of degrees of freedom for the first data set.
        real(real64) :: d2
            !! The measure of degrees of freedom for the second data set.
    contains
        procedure, public :: pdf => fd_pdf
        procedure, public :: cdf => fd_cdf
        procedure, public :: mean => fd_mean
        procedure, public :: median => fd_median
        procedure, public :: mode => fd_mode
        procedure, public :: variance => fd_variance
    end type

    ! distributions_f.f90
    interface
        pure module elemental function fd_pdf(this, x) result(rst)
            !! Computes the probability density function.
            !!
            !! The PDF for a F distribution is given as 
            !! $$ f(x) = 
            !! \sqrt{ \frac{ (d_1 x)^{d_1} d_{2}^{d_2} }{ (d_1 x + d_2)^{d_1 + d_2} } } 
            !! \frac{1}{x \beta \left( \frac{d_1}{2}, \frac{d_2}{2} \right) } $$
            class(f_distribution), intent(in) :: this
                !! The f_distribution object.
            real(real64), intent(in) :: x
                !! The value at which to evaluate the function.
            real(real64) :: rst
                !! The value of the function.
        end function

        pure module elemental function fd_cdf(this, x) result(rst)
            !! Computes the cumulative distribution function.
            !!
            !! The CDF for a F distribution is given as 
            !! $$ F(x) = I_{d_1 x/(d_1 x + d_2)} \left( \frac{d_1}{2}, 
            !! \frac{d_2}{2} \right) $$
            class(f_distribution), intent(in) :: this
                !! The f_distribution object.
            real(real64), intent(in) :: x
                !! The value at which to evaluate the function.
            real(real64) :: rst
                !! The value of the function.
        end function

        pure module function fd_mean(this) result(rst)
            !! Computes the mean of the distribution.
            class(f_distribution), intent(in) :: this
                !! The f_distribution object.
            real(real64) :: rst
                !! The mean.
        end function

        pure module function fd_median(this) result(rst)
            !! Computes the median of the distribution.
            class(f_distribution), intent(in) :: this
                !! The f_distribution object.
            real(real64) :: rst
                !! The median.
        end function

        pure module function fd_mode(this) result(rst)
            !! Computes the mode of the distribution.
            class(f_distribution), intent(in) :: this
                !! The f_distribution object.
            real(real64) :: rst
                !! The mode.
        end function
        
        pure module function fd_variance(this) result(rst)
            !! Computes the variance of the distribution.
            class(f_distribution), intent(in) :: this
                !! The f_distribution object.
            real(real64) :: rst
                !! The variance.
        end function
    end interface

    type, extends(distribution) :: chi_squared_distribution
        !! Defines a Chi-squared distribution.
        integer(int32) :: dof
            !! The number of degrees of freedom.
    contains
        procedure, public :: pdf => cs_pdf
        procedure, public :: cdf => cs_cdf
        procedure, public :: mean => cs_mean
        procedure, public :: median => cs_median
        procedure, public :: mode => cs_mode
        procedure, public :: variance => cs_variance
    end type

    ! distributions_chisquared.f90
    interface
        pure module elemental function cs_pdf(this, x) result(rst)
            !! Computes the probability density function.
            !!
            !! The PDF for a Chi-squared distribution is given as 
            !! $$ f(x) = \frac{x^{k/2 - 1} \exp{-x / 2}} {2^{k / 2} 
            !! \Gamma \left( \frac{k}{2} \right)} $$
            class(chi_squared_distribution), intent(in) :: this
                !! The chi_squared_distribution object.
            real(real64), intent(in) :: x
                !! The value at which to evaluate the function.
            real(real64) :: rst
                !! The value of the function.
        end function

        pure module elemental function cs_cdf(this, x) result(rst)
            !! Computes the cumulative distribution function.
            !!
            !! The CDF for a Chi-squared distribution is given as 
            !! $$ F(x) = \frac{ \gamma \left( \frac{k}{2}, \frac{x}{2} \right) }
            !! { \Gamma \left( \frac{k}{2} \right)} $$
            class(chi_squared_distribution), intent(in) :: this
                !! The chi_squared_distribution object.
            real(real64), intent(in) :: x
                !! The value at which to evaluate the function.
            real(real64) :: rst
                !! The value of the function.
        end function

        pure module function cs_mean(this) result(rst)
            !! Computes the mean of the distribution.
            class(chi_squared_distribution), intent(in) :: this
                !! The chi_squared_distribution object.
            real(real64) :: rst
                !! The mean.
        end function

        pure module function cs_median(this) result(rst)
            !! Computes the median of the distribution.
            class(chi_squared_distribution), intent(in) :: this
                !! The chi_squared_distribution object.
            real(real64) :: rst
                !! The median.
        end function

        pure module function cs_mode(this) result(rst)
            !! Computes the mode of the distribution.
            class(chi_squared_distribution), intent(in) :: this
                !! The chi_squared_distribution object.
            real(real64) :: rst
                !! The mode.
        end function

        pure module function cs_variance(this) result(rst)
            !! Computes the variance of the distribution.
            class(chi_squared_distribution), intent(in) :: this
                !! The chi_squared_distribution object.
            real(real64) :: rst
                !! The variance.
        end function
    end interface

    type, extends(distribution) :: binomial_distribution
        !! Defines a binomial distribution.  The binomial distribution describes
        !! the probability p of getting k successes in n independent trials.
        integer(int32) :: n
            !! The number of independent trials.
        real(real64) :: p
            !! The success probability for each trial.  This parameter must
            !! exist on the set [0, 1].
    contains
        procedure, public :: pdf => bd_pdf
        procedure, public :: cdf => bd_cdf
        procedure, public :: mean => bd_mean
        procedure, public :: median => bd_median
        procedure, public :: mode => bd_mode
        procedure, public :: variance => bd_variance
    end type

    ! distributions_binomial.f90
    interface
        pure module elemental function bd_pdf(this, x) result(rst)
            !! Computes the probability mass function.
            !!
            !! See [this](https://en.wikipedia.org/wiki/Binomial_distribution)
            !! for a description of the probability mass function for this
            !! distribution.
            class(binomial_distribution), intent(in) :: this
                !! The binomial_distribution object.
            real(real64), intent(in) :: x
                !! The value at which to evaluate the function.  This parameter
                !! is the number k successes in the n independent trials.  As
                !! such, this parameter must exist on the set [0, n].
            real(real64) :: rst
                !! The value of the function.
        end function

        pure module elemental function bd_cdf(this, x) result(rst)
            !! Computes the cumulative distribution function.
            !!
            !! See [this](https://en.wikipedia.org/wiki/Binomial_distribution)
            !! for a description of the cumulative distribution function for
            !! this distribution.
            class(binomial_distribution), intent(in) :: this
                !! The binomial_distribution object.
            real(real64), intent(in) :: x
                !! The value at which to evaluate the function.  This parameter
                !! is the number k successes in the n independent trials.  As
                !! such, this parameter must exist on the set [0, n].
            real(real64) :: rst
                !! The value of the function.
        end function

        pure module function bd_mean(this) result(rst)
            !! Computes the mean of the distribution.
            class(binomial_distribution), intent(in) :: this
                !! The binomial_distribution object.
            real(real64) :: rst
                !! The mean.
        end function

        pure module function bd_median(this) result(rst)
            !! Computes the median of the distribution.
            class(binomial_distribution), intent(in) :: this
                !! The binomial_distribution object.
            real(real64) :: rst
                !! The median.
        end function

        pure module function bd_mode(this) result(rst)
            !! Computes the mode of the distribution.
            class(binomial_distribution), intent(in) :: this
                !! The binomial_distribution object.
            real(real64) :: rst
                !! The mode.
        end function

        pure module function bd_variance(this) result(rst)
            !! Computes the variance of the distribution.
            class(binomial_distribution), intent(in) :: this
                !! The chi_squared_distribution object.
            real(real64) :: rst
                !! The variance.
        end function
    end interface

! ******************************************************************************
! GENERAL STATISTICS
! ------------------------------------------------------------------------------
    
    type anova_factor
        !! Defines an ANOVA factor result.
        real(real64) :: dof
            !! The number of degrees of freedome.
        real(real64) :: variance
            !! The estimate of variance.
        real(real64) :: sum_of_squares
            !! The sum of the squares.
        real(real64) :: f_statistic
            !! The F-statistic.
        real(real64) :: probability
            !! The variance probability term.
    end type

    type single_factor_anova_table
        !! Defines a single-factor ANOVA results table. 
        type(anova_factor) :: main_factor
            !! The main, or main factor, results.
        type(anova_factor) :: within_factor
            !! The within-treatement (error) results.
        real(real64) :: total_dof
            !! The total number of degrees of freedom.
        real(real64) :: total_sum_of_squares
            !! The total sum of squares. 
        real(real64) :: total_variance
            !! The total variance estimate.
        real(real64) :: overall_mean
            !! The overall mean value.
    end type

    type two_factor_anova_table
        !! Defines a two-factor ANOVA results table.
        type(anova_factor) :: main_factor_1
            !! The first main-factor results.
        type(anova_factor) :: main_factor_2
            !! The second main-factor results.
        type(anova_factor) :: interaction
            !! The interaction effects.
        type(anova_factor) :: within_factor
            !! The within (error) factor results. 
        real(real64) :: total_dof
            !! The total number of degrees of freedom.
        real(real64) :: total_sum_of_squares
            !! The total sum of squares.
        real(real64) :: total_variance
            !! The total variance estimate.
        real(real64) :: overall_mean
            !! The overall mean value.
    end type
    interface mean
        !! Computes the mean of the values in an array.
        module procedure :: mean_real64
        module procedure :: mean_real32
    end interface

    interface variance
        !! Computes the sample variance of the values in an array.
        !!
        !! The variance computed is the sample variance such that the variance.
        !! $$ s^2 = \frac{\Sigma \left( x_{i} - \bar{x} \right)^2}{n - 1} $$
        module procedure :: variance_real64
        module procedure :: variance_real32
    end interface

    interface standard_deviation
        !! Computes the sample standard deviation of the values in an array.
        !! 
        !! The value computed is the sample standard deviation.
        !! $$ s = \sqrt{ \frac{\Sigma \left( x_{i} - \bar{x} \right)^2}{n - 1} } $$
        module procedure :: standard_deviation_real64
        module procedure :: standard_deviation_real32
    end interface

    interface median
        !! Computes the median of the values in an array.
        module procedure :: median_real64
        module procedure :: median_real32
    end interface

    interface r_squared
        !! Computes the R-squared value for a data set.
        !!
        !! The R-squared value is computed by determining the sum of the squares
        !! of the residuals: 
        !! $$ SS_{res} = \Sigma \left( y_i - f_i \right)^2 $$
        !! The total sum of the squares: 
        !! $$ SS_{tot} = \Sigma \left( y_i - \bar{y} \right)^2 $$. 
        !! The R-squared value is then: 
        !! $$ R^2 = 1 - \frac{SS_{res}}{SS_{tot}} $$.
        !!
        !! See Also:
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Coefficient_of_determination)
        module procedure :: r_squared_real64
        module procedure :: r_squared_real32
    end interface

    interface adjusted_r_squared
        !! Computes the adjusted R-squared value for a data set.
        !!
        !! The adjusted R-squared provides a mechanism for tempering the effects
        !! of extra explanatory variables on the traditional R-squared 
        !! calculation.  It is computed by noting the sample size \( n \) and 
        !! the number of variables \( p \).
        !! $$ \bar{R}^2 = 1 - \left( 1 - R^2 \right) \frac{n - 1}{n - p} $$.
        !!
        !! See Also:
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Coefficient_of_determination#Adjusted_R2)
        module procedure :: adjusted_r_squared_real64
        module procedure :: adjusted_r_squared_real32
    end interface

    interface quantile
        !! Computes the specified quantile of a data set using the SAS 
        !! Method 4.
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Quantile)
        module procedure :: quantile_real64
        module procedure :: quantile_real32
    end interface

    interface t_test_equal_variance
        !! Computes the 2-tailed Student's T-Test for two data sets of 
        !! assumed equivalent variances.
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Student%27s_t-test)
        module procedure :: t_test_equal_var_real64
        module procedure :: t_test_equal_var_real32
    end interface

    interface t_test_unequal_variance
        !! Computes the 2-tailed Student's T-Test for two data sets of 
        !! assumed non-equivalent variances.
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Student%27s_t-test)
        module procedure :: t_test_unequal_var_real64
        module procedure :: t_test_unequal_var_real32
    end interface
    
    interface t_test_paired
        !! Computes the 2-tailed Student's T-Test for two paired data sets.
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Student%27s_t-test)
        module procedure :: t_test_paired_real64
        module procedure :: t_test_paired_real32
    end interface

    interface f_test
        !! Computes the F-test and returns the probability (two-tailed) that
        !! the variances of two data sets are not significantly different.
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/F-test)
        module procedure :: f_test_real64
        module procedure :: f_test_real32
    end interface

    interface anova
        !! Performs an analysis of variance (ANOVA) on the supplied data 
        !! set.
        !!
        !! The following example illustrates a single-factor ANOVA on a 
        !! data set.
        !! ```fortran
        !! program example
        !!     use iso_fortran_env
        !!     use fstats
        !!     implicit none
        !!
        !!     ! Local Variables
        !!     character, parameter :: tab = achar(9)
        !!     real(real64) :: x(10, 2)
        !!     type(single_factor_anova_table) :: tbl
        !!
        !!     ! Define the data
        !!     x = reshape( &
        !!         [ &
        !!             3.086d3, 3.082d3, 3.069d3, 3.072d3, 3.045d3, 3.070d3, 3.079d3, &
        !!             3.050d3, 3.062d3, 3.062d3, 3.075d3, 3.061d3, 3.063d3, 3.038d3, &
        !!             3.070d3, 3.062d3, 3.070d3, 3.049d3, 3.042d3, 3.063d3 &
        !!         ], &
        !!         [10, 2] &
        !!     )
        !!
        !!     ! Perform the ANOVA
        !!     tbl = anova(x)
        !!
        !!     ! Print out the table
        !!     print '(A)', "Description" // tab // "DOF" // tab // "Sum of Sq." // &
        !!         tab // "Variance" // tab // "F-Stat" // tab // "P-Value"
        !!     print '(AF2.0AF5.1AF5.1AF5.3AF5.3)', "Main Factor: " // tab, &
        !!         tbl%main_factor%dof, tab, &
        !!         tbl%main_factor%sum_of_squares, tab // tab, &
        !!         tbl%main_factor%variance, tab // tab, &
        !!         tbl%main_factor%f_statistic, tab, &
        !!         tbl%main_factor%probability
        !!
        !!     print '(AF3.0AF6.1AF5.1)', "Within: " // tab, &
        !!         tbl%within_factor%dof, tab, &
        !!         tbl%within_factor%sum_of_squares, tab // tab, &
        !!         tbl%within_factor%variance
        !!
        !!     print '(AF3.0AF6.1AF5.1)', "Total: " // tab // tab, &
        !!         tbl%total_dof, tab, &
        !!         tbl%total_sum_of_squares, tab // tab, &
        !!         tbl%total_variance
        !!
        !!     print '(AF6.1)', "Overall Mean: ", tbl%overall_mean
        !! end program
        !! ```
        !! The above program produces the following output.
        !! ```text
        !! Description     DOF     Sum of Sq.      Variance        F-Stat  P-Value
        !! Main Factor:    1.      352.8           352.8           2.147   0.160
        !! Within:         18.     2958.2          164.3
        !! Total:          19.     3311.0          174.3
        !! Overall Mean: 3063.5
        !! ```
        !!
        !! See Also
        !! 
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Analysis_of_variance)
        !! - [SPC Excel Single Factor ANOVA](https://www.spcforexcel.com/knowledge/root-cause-analysis/single-factor-anova)
        !! - [SPC Excel Gage R&R](https://www.spcforexcel.com/knowledge/measurement-systems-analysis/anova-gage-rr-part-1)
        !! - [SPC Excel Understanding Regression Statistics](https://www.spcforexcel.com/knowledge/root-cause-analysis/understanding-regression-statistics-part-1)
        !! - [NIST - Two Way ANOVA](https://www.itl.nist.gov/div898/handbook/prc/section4/prc427.htm)
        module procedure :: anova_1_factor
        module procedure :: anova_2_factor
        module procedure :: anova_model_fit
    end interface
    
    interface confidence_interval
        !! Computes the confidence interval for the specified distribution.
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Confidence_interval)
        module procedure :: confidence_interval_real64
        module procedure :: confidence_interval_real32
        module procedure :: confidence_interval_real64_array
        module procedure :: confidence_interval_real32_array
    end interface

    interface trimmed_mean
        !! Computes the trimmed mean of a data set.
        module procedure :: trimmed_mean_real64
        module procedure :: trimmed_mean_real32
    end interface

    interface difference
        !! Computes the difference between elements in an array.
        module procedure :: difference_real64
        module procedure :: difference_real32
    end interface

    ! statistics_implementation.f90
    interface
        pure module function mean_real64(x) result(rst)
            !! Computes the mean of the values in an array.
            real(real64), intent(in) :: x(:)
                !! The array of values to analyze.
            real(real64) :: rst
                !! The result.
        end function

        pure module function mean_real32(x) result(rst)
            !! Computes the mean of the values in an array.
            real(real32), intent(in) :: x(:)
                !! The array of values to analyze.
            real(real32) :: rst
                !! The result.
        end function

        pure module function variance_real64(x) result(rst)
            !! Computes the sample variance of the values in an array.
            real(real64), intent(in) :: x(:)
                !! The array of values to analyze.
            real(real64) :: rst
                !! The result.
        end function

        pure module function variance_real32(x) result(rst)
            !! Computes the sample variance of the values in an array.
            real(real32), intent(in) :: x(:)
                !! The array of values to analyze.
            real(real32) :: rst
                !! The result.
        end function

        pure module function standard_deviation_real64(x) result(rst)
            !! Computes the sample standard deviation of the values in an array.
            real(real64), intent(in) :: x(:)
                !! The array of values to analyze.
            real(real64) :: rst
                !! The result.
        end function

        pure module function standard_deviation_real32(x) result(rst)
            !! Computes the sample standard deviation of the values in an array.
            real(real32), intent(in) :: x(:)
                !! The array of values to analyze.
            real(real32) :: rst
                !! The result.
        end function

        module function median_real64(x) result(rst)
            !! Computes the median of the values in an array.
            real(real64), intent(inout) :: x(:)
                !! The array of values to analyze.  On output, the same values,
                !! but sorted into ascending order.
            real(real64) :: rst
                !! The result.
        end function

        module function median_real32(x) result(rst)
            !! Computes the median of the values in an array.
            real(real32), intent(inout) :: x(:)
                !! The array of values to analyze.  On output, the same values,
                !! but sorted into ascending order.
            real(real32) :: rst
                !! The result.
        end function

        module function r_squared_real64(x, xm, err) result(rst)
            !! Computes the R-squared value for a data set.
            real(real64), intent(in) :: x(:)
                !! An N-element array containing the dependent variables from 
                !! the data set.
            real(real64), intent(in) :: xm(:)
                !! An N-element array containing the corresponding modeled 
                !! values.
            class(errors), intent(inout), optional, target :: err
                !! A mechanism for communicating errors and warnings
                !! to the caller.  Possible warning and error codes are as 
                !! follows.
                !! - FS_NO_ERROR: No errors encountered.
                !! - FS_ARRAY_SIZE_ERROR: Occurs if @p x and @p xm are not the 
                !!   same size.
            real(real64) :: rst
                !! The result.
        end function

        module function r_squared_real32(x, xm, err) result(rst)
            !! Computes the R-squared value for a data set.
            real(real32), intent(in) :: x(:)
                !! An N-element array containing the dependent variables from 
                !! the data set.
            real(real32), intent(in) :: xm(:)
                !! An N-element array containing the corresponding modeled 
                !! values.
            class(errors), intent(inout), optional, target :: err
                !! A mechanism for communicating errors and warnings
                !! to the caller.  Possible warning and error codes are as 
                !! follows.
                !! - FS_NO_ERROR: No errors encountered.
                !! - FS_ARRAY_SIZE_ERROR: Occurs if @p x and @p xm are not the 
                !!   same size.
            real(real32) :: rst
                !! The result.
        end function

        module function adjusted_r_squared_real64(p, x, xm, err) result(rst)
            !! Computes the adjusted R-squared value for a data set.
            integer(int32), intent(in) :: p
                !! The number of variables.
            real(real64), intent(in) :: x(:)
                !! An N-element array containing the dependent variables from 
                !! the data set.
            real(real64), intent(in) :: xm(:)
                !! An N-element array containing the corresponding modeled 
                !! values.
            class(errors), intent(inout), optional, target :: err
                !! A mechanism for communicating errors and warnings
                !! to the caller.  Possible warning and error codes are as 
                !! follows.
                !! - FS_NO_ERROR: No errors encountered.
                !! - FS_ARRAY_SIZE_ERROR: Occurs if x and xm are not the 
                !!   same size.
            real(real64) :: rst
                !! The result.
        end function

        module function adjusted_r_squared_real32(p, x, xm, err) result(rst)
            !! Computes the adjusted R-squared value for a data set.
            integer(int32), intent(in) :: p
                !! The number of variables.
            real(real32), intent(in) :: x(:)
                !! An N-element array containing the dependent variables from 
                !! the data set.
            real(real32), intent(in) :: xm(:)
                !! An N-element array containing the corresponding modeled 
                !! values.
            class(errors), intent(inout), optional, target :: err
                !! A mechanism for communicating errors and warnings
                !! to the caller.  Possible warning and error codes are as 
                !! follows.
                !! - FS_NO_ERROR: No errors encountered.
                !! - FS_ARRAY_SIZE_ERROR: Occurs if @p x and @p xm are not the 
                !!   same size.
            real(real32) :: rst
                !! The result.
        end function

        pure module function quantile_real64(x, q) result(rst)
            !! Computes the specified quantile of a data set using the SAS 
            !! Method 4.
            real(real64), intent(in) :: x(:)
                !! An N-element array containing the data.
            real(real64), intent(in) :: q
                !! The quantile to compute (e.g. 0.25 computes the 25% 
                !! quantile).
            real(real64) :: rst
                !! The result.
        end function

        pure module function quantile_real32(x, q) result(rst)
            !! Computes the specified quantile of a data set using the SAS 
            !! Method 4.
            real(real32), intent(in) :: x(:)
                !! An N-element array containing the data.
            real(real32), intent(in) :: q
                !! The quantile to compute (e.g. 0.25 computes the 25% 
                !! quantile).
            real(real32) :: rst
                !! The result.
        end function

        pure module function confidence_interval_real64(dist, alpha, s, n) &
            result(rst)
            !! Computes the confidence interval for the specified distribution.
            class(distribution), intent(in) :: dist
                !! The distribution object defining the probability distribution
                !! to establish the confidence level.
            real(real64), intent(in) :: alpha
                !! The probability value of interest.
            real(real64), intent(in) :: s
                !! The sample standard deviation.
            integer(int32), intent(in) :: n
                !! The number of samples in the data set.
            real(real64) :: rst
                !! The result.
        end function

        pure module function confidence_interval_real32(dist, alpha, s, n) &
            result(rst)
            class(distribution), intent(in) :: dist
                !! The distribution object defining the probability distribution
                !! to establish the confidence level.
            real(real32), intent(in) :: alpha
                !! The probability value of interest.
            real(real32), intent(in) :: s
                !! The sample standard deviation.
            integer(int32), intent(in) :: n
                !! The number of samples in the data set.
            real(real32) :: rst
                !! The result.
        end function

        pure module function confidence_interval_real64_array(dist, alpha, x) &
            result(rst)
            class(distribution), intent(in) :: dist
                !! The distribution object defining the probability distribution
                !! to establish the confidence level.
            real(real64), intent(in) :: alpha
                !! The probability value of interest.
            real(real64), intent(in) :: x(:)
                !! An N-element array containing the data to analyze.
            real(real64) :: rst
                !! The result.
        end function

        pure module function confidence_interval_real32_array(dist, alpha, x) &
            result(rst)
            class(distribution), intent(in) :: dist
                !! The distribution object defining the probability distribution
                !! to establish the confidence level.
            real(real32), intent(in) :: alpha
                !! The probability value of interest.
            real(real32), intent(in) :: x(:)
                !! An N-element array containing the data to analyze.
            real(real32) :: rst
                !! The result.
        end function
        
        module subroutine t_test_equal_var_real64(x1, x2, stat, p, dof)
            !! Computes the 2-tailed Student's T-Test for two data sets of 
            !! assumed equivalent variances.
            real(real64), intent(in) :: x1(:)
                !! An N-element array containing the first data set.
            real(real64), intent(in) :: x2(:)
                !! An M-element array containing the second data set.
            real(real64), intent(out) :: stat
                !! The Student-'s T-Test statistic.
            real(real64), intent(out) :: p
                !! The probability value that the two samples are likely to
                !! have come from the same two underlying populations that 
                !! have the same mean.
            real(real64), intent(out) :: dof
                !! The degrees of freedom.
        end subroutine

        module subroutine t_test_equal_var_real32(x1, x2, stat, p, dof)
            !! Computes the 2-tailed Student's T-Test for two data sets of 
            !! assumed equivalent variances.
            real(real32), intent(in) :: x1(:)
                !! An N-element array containing the first data set.
            real(real32), intent(in) :: x2(:)
                !! An M-element array containing the second data set.
            real(real32), intent(out) :: stat
                !! The Student-'s T-Test statistic.
            real(real32), intent(out) :: p
                !! The probability value that the two samples are likely to
                !! have come from the same two underlying populations that 
                !! have the same mean.
            real(real32), intent(out) :: dof
                !! The degrees of freedom.
        end subroutine

        module subroutine t_test_unequal_var_real64(x1, x2, stat, p, dof)
            !! Computes the 2-tailed Student's T-Test for two data sets of 
            !! assumed non-equivalent variances.
            real(real64), intent(in) :: x1(:)
                !! An N-element array containing the first data set.
            real(real64), intent(in) :: x2(:)
                !! An M-element array containing the second data set.
            real(real64), intent(out) :: stat
                !! The Student-'s T-Test statistic.
            real(real64), intent(out) :: p
                !! The probability value that the two samples are likely to
                !! have come from the same two underlying populations that 
                !! have the same mean.
            real(real64), intent(out) :: dof
                !! The degrees of freedom.
        end subroutine

        module subroutine t_test_unequal_var_real32(x1, x2, stat, p, dof)
            !! Computes the 2-tailed Student's T-Test for two data sets of 
            !! assumed non-equivalent variances.
            real(real32), intent(in) :: x1(:)
                !! An N-element array containing the first data set.
            real(real32), intent(in) :: x2(:)
                !! An M-element array containing the second data set.
            real(real32), intent(out) :: stat
                !! The Student-'s T-Test statistic.
            real(real32), intent(out) :: p
                !! The probability value that the two samples are likely to
                !! have come from the same two underlying populations that 
                !! have the same mean.
            real(real32), intent(out) :: dof
                !! The degrees of freedom.
        end subroutine

        module subroutine t_test_paired_real64(x1, x2, stat, p, dof, err)
            !! Computes the 2-tailed Student's T-Test for two paired data sets.
            real(real64), intent(in) :: x1(:)
                !! An N-element array containing the first data set.
            real(real64), intent(in) :: x2(:)
                !! An N-element array containing the second data set.
            real(real64), intent(out) :: stat
                !! The Student-'s T-Test statistic.
            real(real64), intent(out) :: p
                !! The probability value that the two samples are likely to
                !! have come from the same two underlying populations that 
                !! have the same mean.
            real(real64), intent(out) :: dof
                !! The degrees of freedom.
            class(errors), intent(inout), optional, target :: err
                !! A mechanism for communicating errors and warnings to the 
                !! caller.  Possible warning and error codes are as follows.
                !! - FS_NO_ERROR: No errors encountered.
                !! - FS_ARRAY_SIZE_ERROR: Occurs if x1 and x2 are not the same 
                !!   length.
        end subroutine

        module subroutine t_test_paired_real32(x1, x2, stat, p, dof, err)
            !! Computes the 2-tailed Student's T-Test for two paired data sets.
            real(real32), intent(in) :: x1(:)
                !! An N-element array containing the first data set.
            real(real32), intent(in) :: x2(:)
                !! An N-element array containing the second data set.
            real(real32), intent(out) :: stat
                !! The Student-'s T-Test statistic.
            real(real32), intent(out) :: p
                !! The probability value that the two samples are likely to
                !! have come from the same two underlying populations that 
                !! have the same mean.
            real(real32), intent(out) :: dof
                !! The degrees of freedom.
            class(errors), intent(inout), optional, target :: err
                !! A mechanism for communicating errors and warnings to the 
                !! caller.  Possible warning and error codes are as follows.
                !! - FS_NO_ERROR: No errors encountered.
                !! - FS_ARRAY_SIZE_ERROR: Occurs if x1 and x2 are not the same 
                !!   length.
        end subroutine

        module subroutine f_test_real64(x1, x2, stat, p, dof1, dof2)
            !! Computes the F-test and returns the probability (two-tailed) that
            !! the variances of two data sets are not significantly different.
            real(real64), intent(in) :: x1(:), x2(:)
            real(real64), intent(out) :: stat, p, dof1, dof2
        end subroutine

        module subroutine f_test_real32(x1, x2, stat, p, dof1, dof2)
            !! Computes the F-test and returns the probability (two-tailed) that
            !! the variances of two data sets are not significantly different.
            real(real32), intent(in) :: x1(:), x2(:)
            real(real32), intent(out) :: stat, p, dof1, dof2
        end subroutine
        
        module function anova_1_factor(x) result(rst)
            !! Performs an analysis of variance (ANOVA) on the supplied data 
            !! set.
            real(real64), intent(in) :: x(:,:)
                !! An M-by-N matrix containing the M replications of the N test 
                !! points of interest.
            type(single_factor_anova_table) :: rst
                !! A single_factor_anova_table instance containing the ANOVA
                !! results.
        end function

        module function anova_2_factor(x) result(rst)
            !! Performs an analysis of variance (ANOVA) on the supplied data 
            !! set.
            real(real64), intent(in) :: x(:,:,:)
                !! An M-by-N-by-K array containing the M replications of the
                !! N first factor results, and the K second factor results.
            type(two_factor_anova_table) :: rst
                !! A @ref two_factor_anova_table instance containing the ANOVA
                !! results.
        end function

        module function anova_model_fit(nmodelparams, ymeas, ymod, err) &
            result(rst)
            !! Performs an analysis of variance (ANOVA) on the supplied data 
            !! set.
            integer(int32), intent(in) :: nmodelparams
                !! The number of model parameters.
            real(real64), intent(in) :: ymeas(:)
                !! An N-element array containing the measured dependent
                !! variable data.
            real(real64), intent(in) :: ymod(:)
                !! An N-element array containing the modeled dependent 
                !! variable data.
            class(errors), intent(inout), optional, target :: err
                !! A mechanism for communicating errors and warnings to the 
                !! caller.  Possible warning and error codes are as follows.
                !! - FS_NO_ERROR: No errors encountered.
                !! - FS_ARRAY_SIZE_ERROR: Occurs if ymeas and ymod are not the 
                !!   same length.
                !! - FS_MEMORY_ERROR: Occurs if a memory error is encountered.
            type(single_factor_anova_table) :: rst
        end function

        module function trimmed_mean_real64(x, p) result(rst)
            !! Computes the trimmed mean of a data set.
            real(real64), intent(inout), dimension(:) :: x
                !! An N-element array containing the data.  On output, the
                !! array is sorted into ascending order.
            real(real64), intent(in), optional :: p
                !! An optional parameter specifying the percentage of values
                !! from either end of the distribution to remove.  The default
                !! is 0.05 such that the bottom 5% and top 5% are removed.
            real(real64) :: rst
                !! The trimmed mean.
        end function

        module function trimmed_mean_real32(x, p) result(rst)
            !! Computes the trimmed mean of a data set.
            real(real32), intent(inout), dimension(:) :: x
                !! An N-element array containing the data.  On output, the
                !! array is sorted into ascending order.
            real(real32), intent(in), optional :: p
                !! An optional parameter specifying the percentage of values
                !! from either end of the distribution to remove.  The default
                !! is 0.05 such that the bottom 5% and top 5% are removed.
            real(real32) :: rst
                !! The trimmed mean.
        end function

        pure module function difference_real64(x) result(rst)
            !! Computes the difference between elements in an array.
            real(real64), intent(in), dimension(:) :: x
                !! The N-element array on which to operate.
            real(real64), allocatable, dimension(:) :: rst
                !! The (N-1)-element array containing the differences between adjacent
                !! elements.
        end function

        pure module function difference_real32(x) result(rst)
            !! Computes the difference between elements in an array.
            real(real32), intent(in), dimension(:) :: x
                !! The N-element array on which to operate.
            real(real32), allocatable, dimension(:) :: rst
                !! The (N-1)-element array containing the differences between adjacent
                !! elements.
        end function
    end interface

! ******************************************************************************
! SPECIAL FUNCTIONS
! ------------------------------------------------------------------------------
    interface beta
        !! Computes the beta function.
        !!
        !! The beta function is related to the gamma function
        !! by the following relationship.
        !! $$ \beta(a,b) = \frac{\Gamma(a) \Gamma(b)}{\Gamma(a + b)} $$.
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Beta_function)
        module procedure :: beta_real64
        module procedure :: beta_real32
    end interface

    interface regularized_beta
        !! Computes the regularized beta function.
        !!
        !! The regularized beta function is defined as the ratio between
        !! the incomplete beta function and the beta function.
        !! $$ I_{x}(a,b) = \frac{\beta(x;a,b)}{\beta(a,b)} $$.
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Beta_function)
        module procedure :: regularized_beta_real64
        module procedure :: regularized_beta_real32
    end interface

    interface incomplete_beta
        !! Computes the incomplete beta function.
        !!
        !! The incomplete beta function is defind as:
        !! $$ \beta(x;a,b) = \int_{0}^{x} t^{a-1} (1 - t)^{b-1} dt $$.
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Beta_function#Incomplete_beta_function)
        module procedure :: incomplete_beta_real64
        module procedure :: incomplete_beta_real32
    end interface

    interface digamma
        !! Computes the digamma function.
        !!
        !! The digamma function is defined as:
        !! $$ \psi(x) = 
        !! \frac{d}{dx}\left( \ln \left( \Gamma \left( x \right) \right) 
        !! \right) $$
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Digamma_function)
        module procedure :: digamma_real64
        module procedure :: digamma_real32
    end interface
    
    interface incomplete_gamma_upper
        !! Computes the upper incomplete gamma function.
        !!
        !! The upper incomplete gamma function is defined as:
        !! $$ \Gamma(a, x) = \int_{x}^{\infty} t^{a-1} e^{-t} \,dt $$
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Incomplete_gamma_function)
        module procedure :: incomplete_gamma_upper_real64
        module procedure :: incomplete_gamma_upper_real32
    end interface

    interface incomplete_gamma_lower
        !! Computes the lower incomplete gamma function.
        !!
        !! The lower incomplete gamma function is defined as:
        !! $$ \gamma(a, x) = \int_{0}^{x} t^{a-1} e^{-t} \,dt $$
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Incomplete_gamma_function)
        module procedure :: incomplete_gamma_lower_real64
        module procedure :: incomplete_gamma_lower_real32
    end interface

    ! special_functions_beta.f90
    interface
        pure elemental module function beta_real64(a, b) result(rst)
            !! Computes the beta function.
            real(real64), intent(in) :: a
                !! The first argument of the function.
            real(real64), intent(in) :: b
                !! The second argument of the function.
            real(real64) :: rst
                !! The value of the beta function at \( a \) and \( b \).
        end function

        pure elemental module function beta_real32(a, b) result(rst)
            real(real32), intent(in) :: a
                !! The first argument of the function.
            real(real32), intent(in) :: b
                !! The second argument of the function.
            real(real32) :: rst
                !! The value of the beta function at \( a \) and \( b \).
        end function

        pure elemental module function regularized_beta_real64(a, b, x) result(rst)
            !! Computes the regularized beta function.
            real(real64), intent(in) :: a
                !! The first argument of the function.
            real(real64), intent(in) :: b
                !! The second argument of the function.
            real(real64), intent(in) :: x
                !! The upper limit of the integration.
            real(real64) :: rst
                !! The value of the regularized beta function.
        end function

        pure elemental module function regularized_beta_real32(a, b, x) result(rst)
            !! Computes the regularized beta function.
            real(real32), intent(in) :: a
                !! The first argument of the function.
            real(real32), intent(in) :: b
                !! The second argument of the function.
            real(real32), intent(in) :: x
                !! The upper limit of the integration.
            real(real32) :: rst
                !! The value of the regularized beta function.
        end function

        pure elemental module function incomplete_beta_real64(a, b, x) result(rst)
            !! Computes the incomplete beta function.
            real(real64), intent(in) :: a
                !! The first argument of the function.
            real(real64), intent(in) :: b
                !! The second argument of the function.
            real(real64), intent(in) :: x
                !! The upper limit of the integration.
            real(real64) :: rst
                !! The value of the incomplete beta function.
        end function

        pure elemental module function incomplete_beta_real32(a, b, x) result(rst)
            !! Computes the incomplete beta function.
            real(real32), intent(in) :: a
                !! The first argument of the function.
            real(real32), intent(in) :: b
                !! The second argument of the function.
            real(real32), intent(in) :: x
                !! The upper limit of the integration.
            real(real32) :: rst
                !! The value of the incomplete beta function.
        end function
    end interface

    ! special_functions_digamma.f90
    interface
        pure elemental module function digamma_real64(x) result(rst)
            !! Computes the digamma function.
            real(real64), intent(in) :: x
                !! The value at which to evaluate the function.
            real(real64) :: rst
                !! The function value.
        end function

        pure elemental module function digamma_real32(x) result(rst)
            !! Computes the digamma function.
            real(real32), intent(in) :: x
                !! The value at which to evaluate the function.
            real(real32) :: rst
                !! The function value.
        end function
    end interface

    ! special_functions_gamma.f90
    interface
        pure elemental module function incomplete_gamma_upper_real64(a, x) &
            result(rst)
            !! Computes the upper incomplete gamma function.
            real(real64), intent(in) :: a
                !! The coefficient value.
            real(real64), intent(in) :: x
                !! The value at which to evaluate the function.
            real(real64) :: rst
                !! The function value.
        end function

        pure elemental module function incomplete_gamma_upper_real32(a, x) &
            result(rst)
            !! Computes the upper incomplete gamma function.
            real(real32), intent(in) :: a
                !! The coefficient value.
            real(real32), intent(in) :: x
                !! The value at which to evaluate the function.
            real(real32) :: rst
                !! The function value.
        end function

        pure elemental module function incomplete_gamma_lower_real64(a, x) &
            result(rst)
            !! Computes the lower incomplete gamma function.
            real(real64), intent(in) :: a
                !! The coefficient value.
            real(real64), intent(in) :: x
                !! The value at which to evaluate the function.
            real(real64) :: rst
                !! The function value.
        end function

        pure elemental module function incomplete_gamma_lower_real32(a, x) &
            result(rst)
            !! Computes the lower incomplete gamma function.
            real(real32), intent(in) :: a
                !! The coefficient value.
            real(real32), intent(in) :: x
                !! The value at which to evaluate the function.
            real(real32) :: rst
                !! The function value.
        end function
    end interface

! ******************************************************************************
! REGRESSION
! ------------------------------------------------------------------------------
    type regression_statistics
       !! A container for regression-related statistical information. 
        real(real64) :: standard_error
            !! The standard error for the model coefficient.
            !!
            !! $$ E_{s}(\beta_{i}) = \sqrt{\sigma^{2} C_{ii}} $$
        real(real64) :: t_statistic
            !! The T-statistic for the model coefficient.
            !!
            !! $$ t_o = \frac{ \beta_{i} }{E_{s}(\beta_{i})} $$
        real(real64) :: probability
            !! The probability that the coefficient is not statistically 
            !! important.  A statistically important coefficient will have a 
            !! low probability (p-value), typically 0.05 or lower; however, a 
            !! p-value of up to ~0.2 may be acceptable dependent upon the 
            !! problem.  Typically any p-value larger than ~0.2 indicates the 
            !! parameter is not statistically important for the model.
            !!
            !! $$ p = t_{|t_o|, df_{residual}} $$
        real(real64) :: confidence_interval
            !! The confidence interval for the parameter at the level 
            !! determined by the regression process.
            !!
            !! $$ c = t_{\alpha, df} E_{s}(\beta_{i}) $$
    end type

    interface coefficient_matrix
        !! Computes the coefficient matrix \( X \) to the linear 
        !! least-squares regression problem of \( X \beta = y \), where 
        !! \( X \) is the coefficient matrix computed here, \( \beta \) is 
        !! the vector of coefficients to be determined, and \( y \) is the 
        !! vector of measured dependent variables.
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Linear_regression)
        module procedure :: coefficient_matrix_real64
        module procedure :: coefficient_matrix_real32
    end interface

    interface covariance_matrix
        !! Computes the covariance matrix \( C \) where 
        !! \( C = \left( X^{T} X \right)^{-1} \) and \( X \) is computed
        !! by coefficient_matrix.
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Covariance_matrix)
        !! - [Wikipedia - Regression](https://en.wikipedia.org/wiki/Linear_regression)
        module procedure :: covariance_matrix_real64
        module procedure :: covariance_matrix_real32
    end interface

    interface linear_least_squares
        !! Computes a linear least-squares regression to fit a set of data.
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Linear_regression)
        !! - [SPC Excel Understanding Regression Statistics](https://www.spcforexcel.com/knowledge/root-cause-analysis/understanding-regression-statistics-part-1)
        module procedure :: linear_least_squares_real64
        module procedure :: linear_least_squares_real32
    end interface

    interface calculate_regression_statistics
        !! Computes statistics for the quality of fit for a regression model.
        module procedure :: calculate_regression_stats_r64
        module procedure :: calculate_regression_stats_r32
    end interface

    ! regression_implementation.f90
    interface
        module subroutine coefficient_matrix_real64(order, intercept, x, c, err)
            !! Computes the coefficient matrix to the linear-least squares
            !! regression problem.
            integer(int32), intent(in) :: order
                !! The order of the equation to fit.  This value must be
                !! at least one (linear equation), but can be higher as desired.
            logical, intent(in) :: intercept
                !! Set to true if the intercept is being computed
                !! as part of the regression; else, false.
            real(real64), intent(in) :: x(:)
                !! An N-element array containing the independent variable
                !! measurement points.
            real(real64), intent(out) :: c(:,:)
                !! An N-by-K matrix where the results will be written.  K
                !! must equal order + 1 in the event intercept is true; 
                !! however, if intercept is false, K must equal order.
            class(errors), intent(inout), optional, target :: err
                !! A mechanism for communicating errors and warnings to the 
                !! caller.  Possible warning and error codes are as follows.
                !! - FS_NO_ERROR: No errors encountered.
                !! - FS_ARRAY_SIZE_ERROR: Occurs if c is not properly sized.
                !! - FS_INVALID_INPUT_ERROR: Occurs if order is less than 1.
        end subroutine

        module subroutine coefficient_matrix_real32(order, intercept, x, c, err)
            !! Computes the coefficient matrix to the linear-least squares
            !! regression problem.
            integer(int32), intent(in) :: order
                !! The order of the equation to fit.  This value must be
                !! at least one (linear equation), but can be higher as desired.
            logical, intent(in) :: intercept
                !! Set to true if the intercept is being computed
                !! as part of the regression; else, false.
            real(real32), intent(in) :: x(:)
                !! An N-element array containing the independent variable
                !! measurement points.
            real(real32), intent(out) :: c(:,:)
                !! An N-by-K matrix where the results will be written.  K
                !! must equal order + 1 in the event intercept is true; 
                !! however, if intercept is false, K must equal order.
            class(errors), intent(inout), optional, target :: err
                !! A mechanism for communicating errors and warnings to the 
                !! caller.  Possible warning and error codes are as follows.
                !! - FS_NO_ERROR: No errors encountered.
                !! - FS_ARRAY_SIZE_ERROR: Occurs if c is not properly sized.
                !! - FS_INVALID_INPUT_ERROR: Occurs if order is less than 1.
        end subroutine

        module subroutine covariance_matrix_real64(x, c, err)
            !! Computes the covariance matrix.
            real(real64), intent(in) :: x(:,:)
                !! An M-by-N matrix containing the formatted independent data
                !!  matrix \( X \) as computed by coefficient_matrix.
            real(real64), intent(out) :: c(:,:)
                !! The N-by-N covariance matrix.
            class(errors), intent(inout), optional, target :: err
                !! A mechanism for communicating errors and warnings to the 
                !! caller.  Possible warning and error codes are as follows.
                !! - FS_NO_ERROR: No errors encountered.
                !! - FS_ARRAY_SIZE_ERROR: Occurs if any of the matrices are not 
                !!      sized correctly.
                !! - FS_MEMORY_ERROR: Occurs if there is a memory allocation 
                !!      error.
        end subroutine

        module subroutine covariance_matrix_real32(x, c, err)
            !! Computes the covariance matrix.
            real(real32), intent(in) :: x(:,:)
                !! An M-by-N matrix containing the formatted independent data
                !!  matrix \( X \) as computed by coefficient_matrix.
            real(real32), intent(out) :: c(:,:)
                !! The N-by-N covariance matrix.
            class(errors), intent(inout), optional, target :: err
                !! A mechanism for communicating errors and warnings to the 
                !! caller.  Possible warning and error codes are as follows.
                !! - FS_NO_ERROR: No errors encountered.
                !! - FS_ARRAY_SIZE_ERROR: Occurs if any of the matrices are not 
                !!      sized correctly.
                !! - FS_MEMORY_ERROR: Occurs if there is a memory allocation 
                !!      error.
        end subroutine

        module subroutine linear_least_squares_real64(order, intercept, x, y, &
            coeffs, ymod, resid, stats, alpha, err)
            !! Computes a linear least-squares regression to fit a set of data.
            integer(int32), intent(in) :: order
                !! The order of the equation to fit.  This value must be at 
                !! least one (linear equation), but can be higher as desired, 
                !! as long as there is sufficient data.
            logical, intent(in) :: intercept
                !! Set to true if the intercept is being computed as part of 
                !! the regression; else, false.
            real(real64), intent(in) :: x(:)
                !! An N-element array containing the independent variable
                !! measurement points.
            real(real64), intent(in) :: y(:)
                !! An N-element array containing the dependent variable
                !! measurement points.
            real(real64), intent(out) :: coeffs(:)
                !! An ORDER+1 element array where the coefficients will
                !! be written.
            real(real64), intent(out) :: ymod(:)
                !! An N-element array where the modeled data will be written.
            real(real64), intent(out) :: resid(:)
                !! An N-element array where the residual error data will be 
                !! written (modeled - actual).
            type(regression_statistics), intent(out), optional :: stats(:)
                !! An M-element array of regression_statistics items where 
                !! M = ORDER + 1 when intercept is set to true; however, if 
                !! intercept is set to false, M = ORDER.
            real(real64), intent(in), optional :: alpha
                !! The significance level at which to evaluate the confidence 
                !! intervals.  The default value is 0.05 such that a 95% 
                !! confidence interval is calculated.
            class(errors), intent(inout), optional, target :: err
                !! A mechanism for communicating errors and warnings to the 
                !! caller.  Possible warning and error codes are as follows.
                !! - FS_NO_ERROR: No errors encountered.
                !! - FS_ARRAY_SIZE_ERROR: Occurs if any of the arrays are not 
                !!      approriately sized.
                !! - FS_INVALID_INPUT_ERROR: Occurs if order is less than 1.
                !! - FS_MEMORY_ERROR: Occurs if there is a memory allocation 
                !!      error.
        end subroutine

        module subroutine linear_least_squares_real32(order, intercept, x, y, &
            coeffs, ymod, resid, stats, alpha, err)
            !! Computes a linear least-squares regression to fit a set of data.
            integer(int32), intent(in) :: order
                !! The order of the equation to fit.  This value must be at 
                !! least one (linear equation), but can be higher as desired, 
                !! as long as there is sufficient data.
            logical, intent(in) :: intercept
                !! Set to true if the intercept is being computed as part of 
                !! the regression; else, false.
            real(real32), intent(in) :: x(:)
                !! An N-element array containing the independent variable
                !! measurement points.
            real(real32), intent(in) :: y(:)
                !! An N-element array containing the dependent variable
                !! measurement points.
            real(real32), intent(out) :: coeffs(:)
                !! An ORDER+1 element array where the coefficients will
                !! be written.
            real(real32), intent(out) :: ymod(:)
                !! An N-element array where the modeled data will be written.
            real(real32), intent(out) :: resid(:)
                !! An N-element array where the residual error data will be 
                !! written (modeled - actual).
            type(regression_statistics), intent(out), optional :: stats(:)
                !! An M-element array of regression_statistics items where 
                !! M = ORDER + 1 when intercept is set to true; however, if 
                !! intercept is set to false, M = ORDER.
            real(real32), intent(in), optional :: alpha
                !! The significance level at which to evaluate the confidence 
                !! intervals.  The default value is 0.05 such that a 95% 
                !! confidence interval is calculated.
            class(errors), intent(inout), optional, target :: err
                !! A mechanism for communicating errors and warnings to the 
                !! caller.  Possible warning and error codes are as follows.
                !! - FS_NO_ERROR: No errors encountered.
                !! - FS_ARRAY_SIZE_ERROR: Occurs if any of the arrays are not 
                !!      approriately sized.
                !! - FS_INVALID_INPUT_ERROR: Occurs if order is less than 1.
                !! - FS_MEMORY_ERROR: Occurs if there is a memory allocation 
                !!      error.
        end subroutine

        module function calculate_regression_stats_r64(resid, params, c, &
            alpha, err) result(rst)
            !! Computes statistics for the quality of fit for a regression 
            !! model.
            real(real64), intent(in) :: resid(:)
                !! An M-element array containing the model residual errors.
            real(real64), intent(in) :: params(:)
                !! An N-element array containing the model parameters.
            real(real64), intent(in) :: c(:,:)
                !! The N-by-N covariance matrix.
            real(real64), intent(in), optional :: alpha
                !! The significance level at which to evaluate the confidence 
                !! intervals.  The default value is 0.05 such that a 95% 
                !! confidence interval is calculated.
            class(errors), intent(inout), optional, target :: err
                !! A mechanism for communicating errors and warnings to the 
                !! caller.  Possible warning and error codes are as follows.
                !! - FS_NO_ERROR: No errors encountered.
                !! - FS_ARRAY_SIZE_ERROR: Occurs if c is not sized correctly.
                !! - FS_INVALID_INPUT_ERROR: Occurs if order is less than 1.
                !! - FS_MEMORY_ERROR: Occurs if there is a memory allocation 
                !!      error.
            type(regression_statistics), allocatable :: rst(:)
                !! A regression_statistics object containing the analysis
                !! results.
        end function

        module function calculate_regression_stats_r32(resid, params, c, &
            alpha, err) result(rst)
            !! Computes statistics for the quality of fit for a regression 
            !! model.
            real(real32), intent(in) :: resid(:)
                !! An M-element array containing the model residual errors.
            real(real32), intent(in) :: params(:)
                !! An N-element array containing the model parameters.
            real(real32), intent(in) :: c(:,:)
                !! The N-by-N covariance matrix.
            real(real32), intent(in), optional :: alpha
                !! The significance level at which to evaluate the confidence 
                !! intervals.  The default value is 0.05 such that a 95% 
                !! confidence interval is calculated.
            class(errors), intent(inout), optional, target :: err
                !! A mechanism for communicating errors and warnings to the 
                !! caller.  Possible warning and error codes are as follows.
                !! - FS_NO_ERROR: No errors encountered.
                !! - FS_ARRAY_SIZE_ERROR: Occurs if c is not sized correctly.
                !! - FS_INVALID_INPUT_ERROR: Occurs if order is less than 1.
                !! - FS_MEMORY_ERROR: Occurs if there is a memory allocation 
                !!      error.
            type(regression_statistics), allocatable :: rst(:)
                !! A regression_statistics object containing the analysis
                !! results.
        end function
    end interface

! ******************************************************************************
! EXPERIMENTAL DESIGN
! ------------------------------------------------------------------------------
    interface get_full_factorial_matrix_size
        !! Computes the appropriate size for a full-factorial design table.
        module procedure :: get_full_factorial_matrix_size_int32
    end interface

    interface full_factorial
        !! Computes a table with values scaled from 1 to N describing a 
        !! full-factorial design.
        !!
        !! ```fortran
        !! program example
        !!     use iso_fortran_env
        !!     use fstats
        !!     implicit none
        !!
        !!     ! Local Variables
        !!     integer(int32) :: i, vars(3), tbl(24, 3)
        !!
        !!     ! Define the number of design points for each of the 3 factors to study
        !!     vars = [2, 4, 3]
        !!
        !!     ! Determine the design table
        !!     call full_factorial(vars, tbl)
        !!
        !!     ! Display the table
        !!     do i = 1, size(tbl, 1)
        !!         print *, tbl(i,:)
        !!     end do
        !! end program
        !! ```
        !! The above program produces the following output.
        !! ```text
        !! 1           1           1
        !! 1           1           2
        !! 1           1           3
        !! 1           2           1
        !! 1           2           2
        !! 1           2           3
        !! 1           3           1
        !! 1           3           2
        !! 1           3           3
        !! 1           4           1
        !! 1           4           2
        !! 1           4           3
        !! 2           1           1
        !! 2           1           2
        !! 2           1           3
        !! 2           2           1
        !! 2           2           2
        !! 2           2           3
        !! 2           3           1
        !! 2           3           2
        !! 2           3           3
        !! 2           4           1
        !! 2           4           2
        !! 2           4           3
        !! ```
        module procedure :: full_factorial_int32
    end interface

    ! experimental_design_implementation.f90
    interface
        module subroutine get_full_factorial_matrix_size_int32(vars, m, n, err)
            !! Computes the appropriate size for a full-factorial design table.
            integer(int32), intent(in) :: vars(:)
                !! An M-element array containing the M factors to study.  Each 
                !! of the M entries to the array is expected to contain the 
                !! number of options for that particular factor to explore.  
                !! This value must be greater than or equal to 1.
            integer(int32), intent(out) :: m
                !! The number of rows for the table.
            integer(int32), intent(out) :: n
                !! The number of columns for the table.
            class(errors), intent(inout), optional, target :: err
                !! A mechanism for communicating errors and warnings to the 
                !! caller.  Possible warning and error codes are as follows.
                !! - FS_NO_ERROR: No errors encountered.
                !! - FS_INVALID_INPUT_ERROR: Occurs if any items in vars are 
                !!      less than 1.
        end subroutine

        module subroutine full_factorial_int32(vars, tbl, err)
            !! Computes a table with values scaled from 1 to N describing a 
            !! full-factorial design.
            integer(int32), intent(in) :: vars(:)
                !! An M-element array containing the M factors to study.  
                !! Each of the M entries to the array is expected to contain 
                !! the number of options for that particular factor to explore. 
                !! This value must be greater than or equal to 1.
            integer(int32), intent(out) :: tbl(:,:)
                !! A table where the design will be written.  Use 
                !! get_full_factorial_matrix_size to determine the appropriate 
                !! table size.
            class(errors), intent(inout), optional, target :: err
                !! A mechanism for communicating errors and warnings to the 
                !! caller.  Possible warning and error codes are as follows.
                !! - FS_NO_ERROR: No errors encountered.
                !! - FS_INVALID_INPUT_ERROR: Occurs if any items in vars are 
                !!      less than 1.
                !! - FS_ARRAY_SIZE_ERROR: Occurs if tbl is not properly sized.
        end subroutine
    end interface

! ******************************************************************************
! NONLINEAR REGRESSION
! ------------------------------------------------------------------------------
    type iteration_controls
        !! Provides a collection of iteration control parameters.
        integer(int32) :: max_iteration_count
            !! Defines the maximum number of iterations allowed.
        integer(int32) :: max_function_evaluations
            !! Defines the maximum number of function evaluations allowed.
        real(real64) :: gradient_tolerance
            !! Defines a tolerance on the gradient of the fitted function.
        real(real64) :: change_in_solution_tolerance
            !! Defines a tolerance on the change in parameter values.
        real(real64) :: residual_tolerance
            !! Defines a tolerance on the metric associated with the residual 
            !! error.
        real(real64) :: iteration_improvement_tolerance
            !! Defines a tolerance to ensure adequate improvement on each 
            !! iteration.
        integer(int32) :: max_iteration_between_updates
            !! Defines how many iterations can pass before a re-evaluation of 
            !! the Jacobian matrix is forced.
    contains
        procedure, public :: set_to_default => lm_set_default_tolerances
    end type

    type convergence_info
        !! Provides information regarding convergence status.
        logical :: converge_on_gradient
            !! True if convergence on the gradient was achieved; else, false.
        real(real64) :: gradient_value
            !! The value of the gradient test parameter.
        logical :: converge_on_solution_change
            !! True if convergence on the change in solution was achieved; else,
            !! false.
        real(real64) :: solution_change_value
            !! The value of the change in solution parameter.
        logical :: converge_on_residual_parameter
            !! True if convergence on the residual error parameter was achieved; 
            !! else, false.
        real(real64) :: residual_value
            !! The value of the residual error parameter.
        logical :: reach_iteration_limit
            !! True if the solution did not converge in the allowed number of 
            !! iterations.
        integer(int32) :: iteration_count
            !! The iteration count.
        logical :: reach_function_evaluation_limit
            !! True if the solution did not converge in the allowed number of
            !! function evaluations.
        integer(int32) :: function_evaluation_count
            !! The function evaluation count.
        logical :: user_requested_stop
            !! True if the user requested the stop; else, false.
    end type

    type lm_solver_options
        !! Options to control the Levenberg-Marquardt solver.
        integer(int32) :: method
            !! The solver method to utilize.
            !! - FS_LEVENBERG_MARQUARDT_UPDATE:
            !! - FS_QUADRATIC_UPDATE:
            !! - FS_NIELSEN_UDPATE:
        real(real64) :: finite_difference_step_size
            !! The step size used for the finite difference calculations of the
            !! Jacobian matrix.
        real(real64) :: damping_increase_factor
            !! The factor to use when increasing the damping parameter.
        real(real64) :: damping_decrease_factor
            !! The factor to use when decreasing the damping parameter.
    contains
        procedure, public :: set_to_default => lm_set_default_settings
    end type

    interface
        subroutine regression_function(xdata, params, resid, stop)
            use iso_fortran_env, only : real64
            real(real64), intent(in), dimension(:) :: xdata, params
            real(real64), intent(out), dimension(:) :: resid
            logical, intent(out) :: stop
        end subroutine

        subroutine iteration_update(iter, funvals, resid, params, step)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: iter
            real(real64), intent(in) :: funvals(:), resid(:), params(:), step(:)
        end subroutine
        
        module subroutine regression_jacobian_1(fun, xdata, params, &
            jac, stop, f0, f1, step, err)
            !! Computes the Jacobian matrix for a nonlinear regression problem.
            procedure(regression_function), intent(in), pointer :: fun
                !! A pointer to the regression_function to evaluate.
            real(real64), intent(in) :: xdata(:)
                !! The M-element array containing x-coordinate data.
            real(real64), intent(in) :: params(:)
                !! The N-element array containing the model parameters.
            real(real64), intent(out) :: jac(:,:)
                !! The M-by-N matrix where the Jacobian will be written.
            logical, intent(out) :: stop
                !! A value that the user can set in fun forcing the
                !! evaluation process to stop prior to completion.
            real(real64), intent(in), optional, target :: f0(:)
                !! An optional M-element array containing the model values
                !!  using the current parameters as defined in m.  This input 
                !! can be used to prevent the routine from performing a 
                !! function evaluation at the model parameter state defined in 
                !! params.
            real(real64), intent(out), optional, target :: f1(:)
                !! An optional M-element workspace array used for function
                !! evaluations.
            real(real64), intent(in), optional :: step
                !! The differentiation step size.  The default is the square 
                !! root of machine precision.
            class(errors), intent(inout), optional, target :: err
                !! A mechanism for communicating errors and warnings to the 
                !! caller.  Possible warning and error codes are as follows.
                !! - FS_NO_ERROR: No errors encountered.
                !! - FS_ARRAY_SIZE_ERROR: Occurs if any of the arrays are not 
                !!      properly sized.
                !! - FS_MEMORY_ERROR: Occurs if there is a memory allocation 
                !!      error.
        end subroutine

        module subroutine nonlinear_least_squares_1(fun, x, y, params, ymod, &
            resid, weights, maxp, minp, stats, alpha, controls, settings, &
            info, status, err)
            !! Performs a nonlinear regression to fit a model using a version
            !! of the Levenberg-Marquardt algorithm.
            procedure(regression_function), intent(in), pointer :: fun
                !! A pointer to the regression_function to evaluate.
            real(real64), intent(in) :: x(:)
                !! The M-element array containing independent data.
            real(real64), intent(in) :: y(:)
                !! The M-element array containing dependent data.
            real(real64), intent(inout) :: params(:)
                !! On input, the N-element array containing the initial estimate
                !! of the model parameters.  On output, the computed model 
                !! parameters.
            real(real64), intent(out) :: ymod(:)
                !! An M-element array where the modeled dependent data will
                !! be written.
            real(real64), intent(out) :: resid(:)
                !! An M-element array where the model residuals will be
                !! written.
            real(real64), intent(in), optional, target :: weights(:)
                !! An optional M-element array allowing the weighting of
                !! individual points.
            real(real64), intent(in), optional, target :: maxp(:)
                !! An optional N-element array that can be used as upper limits 
                !! on the parameter values.  If no upper limit is requested for
                !! a particular parameter, utilize a very large value.  The 
                !! internal default is to utilize huge() as a value.
            real(real64), intent(in), optional, target :: minp(:)
                !! An optional N-element array that can be used as lower limits 
                !! on the parameter values.  If no lower limit is requested for
                !! a particalar parameter, utilize a very large magnitude, but 
                !! negative, value.  The internal default is to utilize -huge() 
                !! as a value.
            type(regression_statistics), intent(out), optional :: stats(:)
                !! An optional N-element array that, if supplied, will be used 
                !! to return statistics about the fit for each parameter.
            real(real64), intent(in), optional :: alpha
                !! The significance level at which to evaluate the confidence 
                !! intervals.  The default value is 0.05 such that a 95% 
                !! confidence interval is calculated.
            type(iteration_controls), intent(in), optional :: controls
                !! An optional input providing custom iteration controls.
            type(lm_solver_options), intent(in), optional :: settings
                !! An optional input providing custom settings for the solver.
            type(convergence_info), intent(out), optional, target :: info
                !! An optional output that can be used to gain information about
                !! the iterative solution and the nature of the convergence.
            procedure(iteration_update), intent(in), pointer, optional :: status
                !! An optional pointer to a routine that can be used to extract
                !! iteration information.
            class(errors), intent(inout), optional, target :: err
                !! A mechanism for communicating errors and warnings to the 
                !! caller.  Possible warning and error codes are as follows.
                !! - FS_NO_ERROR: No errors encountered.
                !! - FS_ARRAY_SIZE_ERROR: Occurs if any of the arrays are not 
                !!      properly sized.
                !! - FS_MEMORY_ERROR: Occurs if there is a memory allocation 
                !!      error.
                !! - FS_UNDERDEFINED_PROBLEM_ERROR: Occurs if the problem posed 
                !!      is underdetetermined (M < N).
                !! - FS_TOLERANCE_TOO_SMALL_ERROR: Occurs if any supplied 
                !!      tolerances are too small to be practical.
                !! - FS_TOO_FEW_ITERATION_ERROR: Occurs if too few iterations 
                !!      are allowed.
        end subroutine

        module subroutine ic_equal(x, y)
            type(iteration_controls), intent(inout) :: x
            type(iteration_controls), intent(in) :: y
        end subroutine

        module subroutine ci_equal(x, y)
            type(convergence_info), intent(inout) :: x
            type(convergence_info), intent(in) :: y
        end subroutine

        module subroutine lso_equal(x, y)
            type(lm_solver_options), intent(inout) :: x
            type(lm_solver_options), intent(in) :: y
        end subroutine

        module subroutine lm_set_default_tolerances(x)
            !! Sets the object to its default values.
            class(iteration_controls), intent(inout) :: x
                !! The iteration_controls object.
        end subroutine

        module subroutine lm_set_default_settings(x)
            !! Sets the object to its default values.
            class(lm_solver_options), intent(inout) :: x
                !! The lm_solver_options object.
        end subroutine
    end interface

    interface jacobian
        !! Computes the Jacobian matrix for a nonlinear regression problem.
        module procedure :: regression_jacobian_1
    end interface

    interface nonlinear_least_squares
        !! Performs a nonlinear regression to fit a model using a version
        !! of the Levenberg-Marquardt algorithm.
        module procedure :: nonlinear_least_squares_1
    end interface

! ******************************************************************************
! ALLAN VARIANCE
! ------------------------------------------------------------------------------

    ! allan.f90
    interface
        module function allan_variance_1(x, dt, err) result(rst)
            !! Computes the Allan variance of a data set.
            real(real64), intent(in), dimension(:) :: x
                !! The N-element data set to analyze.
            real(real64), intent(in), optional :: dt
                !! An optional input specifying the time increment between 
                !! samples in x.  If not specified, this value is set to 1.
            class(errors), intent(inout), optional, target :: err
                !! A mechanism for communicating errors and warnings to the 
                !! caller.  Possible warning and error codes are as follows.
                !! - FS_NO_ERROR: No errors encountered.
                !! - FS_MEMORY_ERROR: Occurs if there is a memory allocation 
                !!      error.
            real(real64), allocatable, dimension(:,:) :: rst
                !! An M-by-2 array containing the results where M is N / 2 - 1
                !! if N is even; else, M is (N - 1) / 2 - 1 if N is odd.  The 
                !! first column contains the averaging times associated with 
                !! the M results stored in the second column.
        end function
    end interface

    interface allan_variance
        !! Computes the Allan variance of a data set.
        !!
        !! Remarks
        !!
        !! This implementation computes the fully overlapped Allan variance 
        !! using the method presented by Yadav et. al.
        !! 
        !! Yadav, Shrikanth & Shastri, Saurav & Chakravarthi, Ghanashyam & Kumar, 
        !! Viraj & Rao, Divya & Agrawal, Vinod. (2018). A Fast, Parallel Algorithm 
        !! for Fully Overlapped Allan Variance and Total Variance for Analysis and 
        !! Modeling of Noise in Inertial Sensors. IEEE Sensors Letters. PP. 1-1. 
        !! 10.1109/LSENS.2018.2829799.
        module procedure :: allan_variance_1
    end interface

! ******************************************************************************
! BOOTSTRAPPING
! ------------------------------------------------------------------------------

    type bootstrap_regression_statistics
        !! A container for regression-related statistical information as 
        !! computed in a bootstrap, or equivalent, calculation.
        real(real64) :: standard_error
            !! The standard error for the model coefficient.
        real(real64) :: t_statistic
            !! The T-statistic for the model coefficient.
            !!
            !! $$ t_o = \frac{ \beta_{i} }{E_{s}(\beta_{i})} $$
        real(real64) :: probability
            !! The probability that the coefficient is not statistically 
            !! important.  A statistically important coefficient will have a 
            !! low probability (p-value), typically 0.05 or lower; however, a 
            !! p-value of up to ~0.2 may be acceptable dependent upon the 
            !! problem.  Typically any p-value larger than ~0.2 indicates the 
            !! parameter is not statistically important for the model.
            !!
            !! $$ p = t_{|t_o|, df_{residual}} $$
        real(real64) :: upper_confidence_interval
            !! The upper limit of the confidence interval for the parameter.
        real(real64) :: lower_confidence_interval
            !! The lower limit of the confidence interval for the parameter.
    end type

    ! bootstrapping.f90
    interface
        module subroutine bs_linear_least_squares_real64(order, intercept, &
            x, y, coeffs, ymod, resid, nsamples, stats, bias, alpha, err)
            integer(int32), intent(in) :: order
            logical, intent(in) :: intercept
            real(real64), intent(in), dimension(:) :: x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:) :: coeffs
            real(real64), intent(out), dimension(:) :: ymod
            real(real64), intent(out), dimension(:) :: resid
            integer(int32), intent(in), optional :: nsamples
            type(bootstrap_regression_statistics), intent(out), optional, &
                dimension(:) :: stats
            real(real64), intent(out), optional, dimension(:) :: bias
            real(real64), intent(in), optional :: alpha
            class(errors), intent(inout), optional, target :: err
        end subroutine
    end interface

    interface bootstrap_linear_least_squares
        module procedure :: bs_linear_least_squares_real64
    end interface

! ******************************************************************************
! SAMPLING
! ------------------------------------------------------------------------------
    ! sampling.f90
    interface
        module function box_muller_sample_real64(mu, sigma) result(rst)
            !! Generates a pair of independent, standard, normally distributed
            !! random values using the Box-Muller transform.
            real(real64), intent(in) :: mu
                !! The mean of the distribution.
            real(real64), intent(in) :: sigma
                !! The standard deviation of the distribution.
            real(real64) :: rst(2)
                !! The pair of random values.
        end function

        module function box_muller_array_real64(mu, sigma, n) result(rst)
            !! Generates an array of normally distributed random values sampled
            !! by the Box-Muller transform.
            real(real64), intent(in) :: mu
                !! The mean of the distribution.
            real(real64), intent(in) :: sigma
                !! The standard deviation of the distribution.
            integer(int32), intent(in) :: n
                !! The number of Box-Muller pairs to generate.
            real(real64), allocatable, dimension(:) :: rst
                !! A 2N-element array containing the N Box-Muller pairs.
        end function

        module function rejection_sample(tdist, n, xmin, xmax) result(rst)
            !! Uses rejection sampling to randomly sample a target distribution.
            class(distribution), intent(in) :: tdist
                !! The distribution to sample
            integer(int32), intent(in) :: n
                !! The number of samples to make.
            real(real64), intent(in) :: xmin
                !! The minimum range to explore.
            real(real64), intent(in) :: xmax
                !! The maximum range to explore.
            real(real64), allocatable, dimension(:) :: rst
                !! An N-element array containing the N samples from the 
                !! distribution.
        end function
    end interface

    interface box_muller_sample
        !! Generates random, normally distributed values via the Box-Muller 
        !! transform.
        module procedure :: box_muller_sample_real64
        module procedure :: box_muller_array_real64
    end interface
end module
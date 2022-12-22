! fstats.f90

module fstats
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
    public :: FS_NO_ERROR
    public :: FS_ARRAY_SIZE_ERROR
    public :: FS_INVALID_INPUT_ERROR
    public :: FS_OUT_OF_MEMORY_ERROR

! ******************************************************************************
! ERROR CODES
! ------------------------------------------------------------------------------
    integer(int32), parameter :: FS_NO_ERROR = 0
    integer(int32), parameter :: FS_ARRAY_SIZE_ERROR = 10000
    integer(int32), parameter :: FS_INVALID_INPUT_ERROR = 10001
    integer(int32), parameter :: FS_OUT_OF_MEMORY_ERROR = 10002

! ******************************************************************************
! DISTRIBUTIONS
! ------------------------------------------------------------------------------
    !> Defines a probability distribution.
    type, abstract :: distribution
    contains
        !> Computes the probability density function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function pdf(class(distribution) this, real(real64 )x)
        !! @endcode
        !! 
        !! @param[in] this The distribution object.
        !! @param[in] x The value at which to evaluate the function.
        !! @return The value of the function.
        procedure(distribution_function), deferred, pass :: pdf
        !> Computes the cumulative distribution function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function cdf(class(distribution) this, real(real64 )x)
        !! @endcode
        !! 
        !! @param[in] this The distribution object.
        !! @param[in] x The value at which to evaluate the function.
        !! @return The value of the function.
        procedure(distribution_function), deferred, pass :: cdf
        !> Computes the mean of the distribution.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function mean(class(distribution) this)
        !! @endcode
        !!
        !! @param[in] this The distribution object.
        !! @return The mean value.
        procedure(distribution_property), deferred, pass :: mean
        !> Computes the median of the distribution.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function median(class(distribution) this)
        !! @endcode
        !!
        !! @param[in] this The distribution object.
        !! @return The median value.
        procedure(distribution_property), deferred, pass :: median
        !> Computes the mode of the distribution.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function mode(class(distribution) this)
        !! @endcode
        !!
        !! @param[in] this The distribution object.
        !! @return The mode value.
        procedure(distribution_property), deferred, pass :: mode
        !> Computes the variance of the distribution.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function variance(class(distribution) this)
        !! @endcode
        !!
        !! @param[in] this The distribution object.
        !! @return The variance value.
        procedure(distribution_property), deferred, pass :: variance
    end type

    interface
        pure elemental function distribution_function(this, x) result(rst)
            use iso_fortran_env, only : real64
            import distribution
            class(distribution), intent(in) :: this
            real(real64), intent(in) :: x
            real(real64) :: rst
        end function

        pure function distribution_property(this) result(rst)
            use iso_fortran_env, only : real64
            import distribution
            class(distribution), intent(in) :: this
            real(real64) :: rst
        end function
    end interface

    !> Defines Student's T-Distribution.
    type, extends(distribution) :: t_distribution
        !> The number of degrees of freedom.
        real(real64) :: dof
    contains
        !> Computes the probability density function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function pdf(class(t_distribution) this, real(real64 )x)
        !! @endcode
        !! 
        !! @param[in] this The t_distribution object.
        !! @param[in] x The value at which to evaluate the function.
        !! @return The value of the function.
        !!
        !! @remarks The PDF for Student's T-Distribution is given as 
        !! \f$ f(t) = \frac{ \Gamma \left( \frac{\nu + 1}{2} \right) }{ \sqrt{\nu \pi} \Gamma \left( \frac{\nu}{2} \right) } \left( 1 + \frac{t^2}{\nu} \right)^{-(\nu + 1) / 2} \f$
        !! where \f$ \nu \f$ is the number of degrees of freedom.
        procedure, public :: pdf => td_pdf
        !> Computes the cumulative distribution function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function cdf(class(t_distribution) this, real(real64 )x)
        !! @endcode
        !! 
        !! @param[in] this The t_distribution object.
        !! @param[in] x The value at which to evaluate the function.
        !! @return The value of the function.
        !!
        !! @remarks The CDF for Student's T-Distribution is given as
        !! \f$ F(t) = \int_{-\infty}^{t} f(u) \,du = 1 - \frac{1}{2} I_{x(t)} \left( \frac{\nu}{2}, \frac{1}{2} \right) \f$
        !! where \f$ x(t) = \frac{\nu}{\nu + t^2} \f$,
        !! \f$ I \f$ is the 
        !! <a href = "https://en.wikipedia.org/wiki/Beta_function#Incomplete_beta_function">
        !! regularized incomplete beta function</a>, and \f$ \nu \f$ is the 
        !! number of degrees of freedom.  This formula is valid for 
        !! \f$ t > 0 \f$.  For values of \f$ t < 0 \f$, the symmetry of the 
        !! function can be exploited.
        procedure, public :: cdf => td_cdf
        !> Computes the mean of the distribution.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function mean(class(t_distribution) this)
        !! @endcode
        !!
        !! @param[in] this The t_distribution object.
        !! @return The mean value.
        procedure, public :: mean => td_mean
        !> Computes the median of the distribution.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function median(class(t_distribution) this)
        !! @endcode
        !!
        !! @param[in] this The t_distribution object.
        !! @return The median value.
        procedure, public :: median => td_median
        !> Computes the mode of the distribution.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function mode(class(t_distribution) this)
        !! @endcode
        !!
        !! @param[in] this The t_distribution object.
        !! @return The mode value.
        procedure, public :: mode => td_mode
        !> Computes the variance of the distribution.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function variance(class(t_distribution) this)
        !! @endcode
        !!
        !! @param[in] this The t_distribution object.
        !! @return The variance value.
        procedure, public :: variance => td_variance
    end type

    ! distributions_t.f90
    interface
        pure module elemental function td_pdf(this, x) result(rst)
            class(t_distribution), intent(in) :: this
            real(real64), intent(in) :: x
            real(real64) :: rst
        end function

        pure module elemental function td_cdf(this, x) result(rst)
            class(t_distribution), intent(in) :: this
            real(real64), intent(in) :: x
            real(real64) :: rst
        end function

        pure module function td_mean(this) result(rst)
            class(t_distribution), intent(in) :: this
            real(real64) :: rst
        end function

        pure module function td_median(this) result(rst)
            class(t_distribution), intent(in) :: this
            real(real64) :: rst
        end function

        pure module function td_mode(this) result(rst)
            class(t_distribution), intent(in) :: this
            real(real64) :: rst
        end function

        pure module function td_variance(this) result(rst)
            class(t_distribution), intent(in) :: this
            real(real64) :: rst
        end function
    end interface

    !> Defines a normal distribution.
    type, extends(distribution) :: normal_distribution
        !> The standard deviation of the distribution.
        real(real64) :: standard_deviation
        !> The mean value of the distribution.
        real(real64) :: mean_value
    contains
        !> Computes the probability density function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function pdf(class(normal_distribution) this, real(real64 )x)
        !! @endcode
        !! 
        !! @param[in] this The normal_distribution object.
        !! @param[in] x The value at which to evaluate the function.
        !! @return The value of the function.
        !!
        !! @remarks The PDF for a normal distribution is given as 
        !! \f$ f(x) = \frac{1}{\sigma \sqrt{2 \pi}} \exp \left(-\frac{1}{2} \left( \frac{x - \mu}{\sigma} \right)^2 \right) \f$.
        procedure, public :: pdf => nd_pdf
        !> Computes the cumulative distribution function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function cdf(class(normal_distribution) this, real(real64 )x)
        !! @endcode
        !! 
        !! @param[in] this The normal_distribution object.
        !! @param[in] x The value at which to evaluate the function.
        !! @return The value of the function.
        !!
        !! @remarks The CDF for a normal distribution is given as 
        !! \f$ F(x) = \frac{1}{2} \left( 1 + erf \left( \frac{x - \mu}{\sigma \sqrt{2}} \right) \right) \f$
        !! where \f$ erf \f$ is the 
        !! <a href="https://en.wikipedia.org/wiki/Error_function">error function</a>.
        procedure, public :: cdf => nd_cdf
        !> Computes the mean of the distribution.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function mean(class(normal_distribution) this)
        !! @endcode
        !!
        !! @param[in] this The normal_distribution object.
        !! @return The mean value.
        procedure, public :: mean => nd_mean
        !> Computes the median of the distribution.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function median(class(normal_distribution) this)
        !! @endcode
        !!
        !! @param[in] this The normal_distribution object.
        !! @return The median value.
        procedure, public :: median => nd_median
        !> Computes the mode of the distribution.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function mode(class(normal_distribution) this)
        !! @endcode
        !!
        !! @param[in] this The normal_distribution object.
        !! @return The mode value.
        procedure, public :: mode => nd_mode
        !> Computes the variance of the distribution.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function variance(class(normal_distribution) this)
        !! @endcode
        !!
        !! @param[in] this The normal_distribution object.
        !! @return The variance value.
        procedure, public :: variance => nd_variance
        !> Standardizes the normal distribution to a mean of 0 and a standard
        !! deviation of 1.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine standardize(class(normal_distribution) this)
        !! @endcode
        !!
        !! @param[in,out] this The normal_distribution object.
        procedure, public :: standardize => nd_standardize
    end type

    ! distributions_normal.f90
    interface
        pure module elemental function nd_pdf(this, x) result(rst)
            class(normal_distribution), intent(in) :: this
            real(real64), intent(in) :: x
            real(real64) :: rst
        end function

        pure module elemental function nd_cdf(this, x) result(rst)
            class(normal_distribution), intent(in) :: this
            real(real64), intent(in) :: x
            real(real64) :: rst
        end function

        pure module function nd_mean(this) result(rst)
            class(normal_distribution), intent(in) :: this
            real(real64) :: rst
        end function

        pure module function nd_median(this) result(rst)
            class(normal_distribution), intent(in) :: this
            real(real64) :: rst
        end function

        pure module function nd_mode(this) result(rst)
            class(normal_distribution), intent(in) :: this
            real(real64) :: rst
        end function

        pure module function nd_variance(this) result(rst)
            class(normal_distribution), intent(in) :: this
            real(real64) :: rst
        end function

        module subroutine nd_standardize(this)
            class(normal_distribution), intent(inout) :: this
        end subroutine
    end interface

    !> Defines an F-distribution.
    type, extends(distribution) :: f_distribution
        !> The measure of degrees of freedom for the first data set.
        real(real64) :: d1
        !> The measure of degrees of freedom for the second data set.
        real(real64) :: d2
    contains
        !> Computes the probability density function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function pdf(class(f_distribution) this, real(real64 )x)
        !! @endcode
        !! 
        !! @param[in] this The f_distribution object.
        !! @param[in] x The value at which to evaluate the function.
        !! @return The value of the function.
        !!
        !! @remarks The PDF for a F distribution is given as 
        !! \f$ f(x) = 
        !! \sqrt{ \frac{ (d_1 x)^{d_1} d_{2}^{d_2} }{ (d_1 x + d_2)^{d_1 + d_2} } } 
        !! \frac{1}{x \beta \left( \frac{d_1}{2}, \frac{d_2}{2} \right) } \f$
        !! where \f$ \beta \f$ is the
        !! <a href="https://en.wikipedia.org/wiki/Beta_function">beta 
        !! function</a>.
        procedure, public :: pdf => fd_pdf
        !> Computes the cumulative distribution function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function cdf(class(f_distribution) this, real(real64 )x)
        !! @endcode
        !! 
        !! @param[in] this The f_distribution object.
        !! @param[in] x The value at which to evaluate the function.
        !! @return The value of the function.
        !!
        !! @remarks The CDF for a F distribution is given as 
        !! \f$ F(x) = I_{d_1 x/(d_1 x + d_2)} \left( \frac{d_1}{2}, \frac{d_2}{2} \right) \f$
        !! where \f$ I \f$ is the 
        !! <a href = "https://en.wikipedia.org/wiki/Beta_function#Incomplete_beta_function">
        !! regularized incomplete beta function</a>.
        procedure, public :: cdf => fd_cdf
        !> Computes the mean of the distribution.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function mean(class(f_distribution) this)
        !! @endcode
        !!
        !! @param[in] this The f_distribution object.
        !! @return The mean value.
        procedure, public :: mean => fd_mean
        !> Computes the median of the distribution.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function median(class(f_distribution) this)
        !! @endcode
        !!
        !! @param[in] this The f_distribution object.
        !! @return The median value.
        procedure, public :: median => fd_median
        !> Computes the mode of the distribution.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function mode(class(f_distribution) this)
        !! @endcode
        !!
        !! @param[in] this The f_distribution object.
        !! @return The mode value.
        procedure, public :: mode => fd_mode
        !> Computes the variance of the distribution.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function variance(class(f_distribution) this)
        !! @endcode
        !!
        !! @param[in] this The f_distribution object.
        !! @return The variance value.
        procedure, public :: variance => fd_variance
    end type

    ! distributions_f.f90
    interface
        pure module elemental function fd_pdf(this, x) result(rst)
            class(f_distribution), intent(in) :: this
            real(real64), intent(in) :: x
            real(real64) :: rst
        end function

        pure module elemental function fd_cdf(this, x) result(rst)
            class(f_distribution), intent(in) :: this
            real(real64), intent(in) :: x
            real(real64) :: rst
        end function

        pure module function fd_mean(this) result(rst)
            class(f_distribution), intent(in) :: this
            real(real64) :: rst
        end function

        pure module function fd_median(this) result(rst)
            class(f_distribution), intent(in) :: this
            real(real64) :: rst
        end function

        pure module function fd_mode(this) result(rst)
            class(f_distribution), intent(in) :: this
            real(real64) :: rst
        end function
        
        pure module function fd_variance(this) result(rst)
            class(f_distribution), intent(in) :: this
            real(real64) :: rst
        end function
    end interface

! ******************************************************************************
! GENERAL STATISTICS
! ------------------------------------------------------------------------------
    !> Defines an ANOVA factor result.
    type anova_factor
        !> The number of degrees of freedome.
        real(real64) :: dof
        !> The estimate of variance.
        real(real64) :: variance
        !> The sum of the squares.
        real(real64) :: sum_of_squares
        !> The F-statistic.
        real(real64) :: f_statistic
        !> The variance probability term.
        real(real64) :: probability
    end type

    !> Defines a single-factor ANOVA results table.
    type single_factor_anova_table
        !> The main, or main factor, results.
        type(anova_factor) :: main_factor
        !> The within-treatement (error) results.
        type(anova_factor) :: within_factor
        !> The total number of degrees of freedom.
        real(real64) :: total_dof
        !> The total sum of squares.
        real(real64) :: total_sum_of_squares
        !> The total variance estimate.
        real(real64) :: total_variance
        !> The overall mean value.
        real(real64) :: overall_mean
    end type

    !> Defines a two-factor ANOVA results table.
    type two_factor_anova_table
        !> The first main-factor results.
        type(anova_factor) :: main_factor_1
        !> The second main-factor results.
        type(anova_factor) :: main_factor_2
        !> The interaction effects.
        type(anova_factor) :: interaction
        !> The within (error) factor results.
        type(anova_factor) :: within_factor
        !> The total number of degrees of freedom.
        real(real64) :: total_dof
        !> The total sum of squares.
        real(real64) :: total_sum_of_squares
        !> The total variance estimate.
        real(real64) :: total_variance
        !> The overall mean value.
        real(real64) :: overall_mean
    end type

    !> Computes the mean of the values in an array.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! real(real64) function mean(real(real64) x(:))
    !! real(real32) function mean(real(real32) x(:))
    !! @endcode
    !!
    !! @param[in] x The N-element array on which to operate.
    !! @return The result.
    interface mean
        module procedure :: mean_real64
        module procedure :: mean_real32
    end interface

    !> Computes the sample variance of the values in an array.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! real(real64) function variance(real(real64) x(:))
    !! real(real32) function variance(real(real32) x(:))
    !! @endcode
    !!
    !! @param[in] x The N-element array on which to operate.
    !! @return The result.
    !!
    !! @par Remarks
    !! The variance computed is the sample variance such that the variance
    !! \f$ s^2 \f$ is computed as \f$ s^2 = \frac{\Sigma \left( x_{i} - \bar{x} \right)^2}{n - 1} \f$.
    interface variance
        module procedure :: variance_real64
        module procedure :: variance_real32
    end interface

    !> Computes the sample standard deviation of the values in an array.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! real(real64) function standard_deviation(real(real64) x(:))
    !! real(real32) function standard_deviation(real(real32) x(:))
    !! @endcode
    !!
    !! @param[in] x The N-element array on which to operate.
    !! @return The result.
    !!
    !! @par Remarks
    !! The variance computed is the sample standard deviation such that the 
    !! standard deviation \f$ s \f$ is computed as 
    !! \f$ s = \sqrt{ \frac{\Sigma \left( x_{i} - \bar{x} \right)^2}{n - 1} } \f$.
    interface standard_deviation
        module procedure :: standard_deviation_real64
        module procedure :: standard_deviation_real32
    end interface

    !> Computes the median of the vlaues in an array.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! real(real64) function median(real(real64) x(:))
    !! real(real32) function median(real(real32) x(:))
    !! @endcode
    !!
    !! @param[in,out] On input, the N-element array on which to operate.  On
    !!  output, the same array but sorted into ascending order.
    !! @return The result.
    interface median
        module procedure :: median_real64
        module procedure :: median_real32
    end interface

    !> Computes the R-squared value for a data set.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! real(real64) function r_squared(real(real64) x(:), real(real64) xm(:), optional class(errors) err)
    !! real(real32) function r_squared(real(real32) x(:), real(real32) xm(:), optional class(errors) err)
    !! @endcode
    !!
    !! @param[in] x An N-element array containing the dependent variables 
    !!  from the data set.
    !! @param[in] xm An N-element array containing the corresponding modeled
    !!  values.
    !! @param[in,out] err A mechanism for communicating errors and warnings
    !!  to the caller.  Possible warning and error codes are as follows.
    !! - FS_NO_ERROR: No errors encountered.
    !! - FS_ARRAY_SIZE_ERROR: Occurs if @p x and @p xm are not the same size.
    !!
    !! @return The result.
    !!
    !! @par Remarks
    !! The R-squared value is computed by determining the sum of the squares
    !! of the residuals: \f$ SS_{res} = \Sigma \left( y_i - f_i \right)^2 \f$,
    !! and the total sum of the squares: 
    !! \f$ SS_{tot} = \Sigma \left( y_i - \bar{y} \right)^2 \f$.  The R-squared
    !! value is then: \f$ R^2 = 1 - \frac{SS_{res}}{SS_{tot}} \f$.
    !!
    !! @par See Also:
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Coefficient_of_determination#Adjusted_R2)
    interface r_squared
        module procedure :: r_squared_real64
        module procedure :: r_squared_real32
    end interface

    !> Computes the adjusted R-squared value for a data set.
    !!
    !! @par Syntax
    !! @code{.f90}
    !!
    !! @endcode
    !!
    !! @param[in] p The number of model parameters.
    !! @param[in] x An N-element array containing the dependent variables 
    !!  from the data set.
    !! @param[in] xm An N-element array containing the corresponding modeled
    !!  values.
    !! @param[in,out] err A mechanism for communicating errors and warnings
    !!  to the caller.  Possible warning and error codes are as follows.
    !! - FS_NO_ERROR: No errors encountered.
    !! - FS_ARRAY_SIZE_ERROR: Occurs if @p x and @p xm are not the same size.
    !!
    !! @return The result.
    !!
    !! @par
    !! The adjusted R-squared provides a mechanism for tempering the effects
    !! of extra explanatory variables on the traditional R-squared calculation.
    !! It is computed by noting the sample size \f$ n \f$ and the number of
    !! variables \f$ p \f$: 
    !! \f$ \bar{R}^2 = 1 - \left( 1 - R^2 \right) \frac{n - 1}{n - p} \f$.
    !!
    !! @par See Also:
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Coefficient_of_determination#Adjusted_R2)
    !! - @ref r_squared
    interface adjusted_r_squared
        module procedure :: adjusted_r_squared_real64
        module procedure :: adjusted_r_squared_real32
    end interface

    !> Computes the specified quantile of a data set using the SAS Method 4.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! real(real64) quantile(real(real64) x, real(real64) q)
    !! real(real32) quantile(real(real32) x, real(real32) q)
    !! @endcode
    !! 
    !! @param[in] x An N-element array containing the data.
    !! @param[in] q The quantile to compute (e.g. 0.25 computes the 25% 
    !!  quantile).
    !!
    !! @return The result.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Quantile)
    interface quantile
        module procedure :: quantile_real64
        module procedure :: quantile_real32
    end interface

    !> Computes the 2-tailed Student's T-Test for two data sets of assumed 
    !! equivalent variances.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! subroutine t_test_equal_variance(real(real64) x1(:), real(real64) x2(:), real(real64) stat, real(real64) p, real(real64) dof)
    !! subroutine t_test_equal_variance(real(real32) x1(:), real(real32) x2(:), real(real32) stat, real(real32) p, real(real32) dof)
    !! @endcode
    !!
    !! @param[in] x1 An N-element array containing the first data set.
    !! @param[in] x2 An M-element array containing the second data set.
    !! @param[out] stat The Student-'s T-Test statistic.
    !! @param[out] p The probability value that the two samples are likely to
    !!  have come from the same two underlying populations that have the same
    !!  mean.
    !! @param[out] dof The degrees of freedom.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Student%27s_t-test)
    interface t_test_equal_variance
        module procedure :: t_test_equal_var_real64
        module procedure :: t_test_equal_var_real32
    end interface
    
    !> Computes the 2-tailed Student's T-Test for two data sets of assumed 
    !! non-equivalent variances.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! subroutine t_test_unequal_variance(real(real64) x1(:), real(real64) x2(:), real(real64) stat, real(real64) p, real(real64) dof)
    !! subroutine t_test_unequal_variance(real(real32) x1(:), real(real32) x2(:), real(real32) stat, real(real32) p, real(real32) dof)
    !! @endcode
    !!
    !! @param[in] x1 An N-element array containing the first data set.
    !! @param[in] x2 An M-element array containing the second data set.
    !! @param[out] stat The Student-'s T-Test statistic.
    !! @param[out] p The probability value that the two samples are likely to
    !!  have come from the same two underlying populations that have the same
    !!  mean.
    !! @param[out] dof The degrees of freedom.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Student%27s_t-test)
    interface t_test_unequal_variance
        module procedure :: t_test_unequal_var_real64
        module procedure :: t_test_unequal_var_real32
    end interface

    !> Computes the 2-tailed Student's T-Test for two paired data sets.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! subroutine t_test_paired(real(real64) x1(:), real(real64) x2(:), real(real64) stat, real(real64) p, real(real64) dof, optional class(errors) err)
    !! subroutine t_test_paired(real(real32) x1(:), real(real32) x2(:), real(real32) stat, real(real32) p, real(real32) dof, optional class(errors) err)
    !! @endcode
    !!
    !! @param[in] x1 An N-element array containing the first data set.
    !! @param[in] x2 An N-element array containing the second data set.
    !! @param[out] stat The Student-'s T-Test statistic.
    !! @param[out] p The probability value that the two samples are likely to
    !!  have come from the same two underlying populations that have the same
    !!  mean.
    !! @param[out] dof The degrees of freedom.
    !! @param[in,out] err A mechanism for communicating errors and warnings
    !!  to the caller.  Possible warning and error codes are as follows.
    !! - FS_NO_ERROR: No errors encountered.
    !! - FS_ARRAY_SIZE_ERROR: Occurs if @p x1 and @p x2 are not the same length.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Student%27s_t-test)
    interface t_test_paired
        module procedure :: t_test_paired_real64
        module procedure :: t_test_paired_real32
    end interface

    !> Computes the F-test and returns the probability (two-tailed) that the
    !! variances of two data sets are not significantly different.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! subroutine f_test(real(real64) x1(:), real(real64) x2(:), real(real64) stat, real(real64) p, real(real64) dof1, real(real64) dof2)
    !! subroutine f_test(real(real32) x1(:), real(real32) x2(:), real(real32) stat, real(real32) p, real(real32) dof1, real(real32) dof2)
    !! @endcode
    !!
    !! @param[in] x1 An N-element array containing the first data set.
    !! @param[in] x2 An M-element array containing the second data set.
    !! @param[out] stat The F-Test statistic.
    !! @param[out] p The probability value that the two data sets have an 
    !!  equivalent variance.
    !! @param[out] dof1 The degrees of freedom of the first data set.
    !! @param[out] dof2 The degrees of freedom of the second data set.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/F-test)
    interface f_test
        module procedure :: f_test_real64
        module procedure :: f_test_real32
    end interface

    !> Performs an analysis of variance (ANOVA) on the supplied data set.
    !!
    !! @par Syntax - Single-Factor ANOVA
    !! @code{.f90}
    !! type(single_factor_anova_table) function anova(real(real64) x(:,:))
    !! @endcode
    !!
    !! @param[in] x An M-by-N matrix containing the M replications of the N test 
    !!  points of interest.
    !! @return A @ref single_factor_anova_table instance containing the ANOVA
    !!  results.
    !!
    !! @par Syntax - Two-Factor ANOVA
    !! @code{.f90}
    !! type(two_factor_anova_table) function anova(real(real64) x(:,:,:))
    !! @endcode
    !!
    !! @param[in] x An M-by-N-by-K array containing the M replications of the
    !!  N first factor results, and the K second factor results.
    !! @return A @ref two_factor_anova_table instance containing the ANOVA
    !!  results.
    !!
    !! @par Syntax - Model Fit ANOVA
    !! @code{.f90}
    !! type(single_factor_anova_table) function anova(integer(int32) nmodelparams, real(real64) ymeas(:), real(real64) ymod(:), optional class(errors) err)
    !! @endcode
    !!
    !! @param[in] nmodelparams: The number of model parameters.
    !! @param[in] ymeas An N-element array containing the measured dependent
    !!  variable data.
    !! @param[in] ymod An N-element array containing the modeled dependent 
    !!  variable data.
    !! @param[in,out] err A mechanism for communicating errors and warnings
    !!  to the caller.  Possible warning and error codes are as follows.
    !! - FS_NO_ERROR: No errors encountered.
    !! - FS_ARRAY_SIZE_ERROR: Occurs if @p ymeas and @p ymod are not the same 
    !!      length.
    !! - FS_OUT_OF_MEMORY_ERROR: Occurs if a memory error is encountered.
    !! @return A @ref single_factor_anova_table instance containing the ANOVA
    !!  results.
    !!
    !! @par Example
    !! The following example illustrates a single-factor ANOVA.
    !! @code{.f90}
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
    !! @endcode
    !! The above program produces the following output.
    !! @code{.txt}
    !! Description     DOF     Sum of Sq.      Variance        F-Stat  P-Value
    !! Main Factor:    1.      352.8           352.8           2.147   0.160
    !! Within:         18.     2958.2          164.3
    !! Total:          19.     3311.0          174.3
    !! Overall Mean: 3063.5
    !! @endcode
    !! For comparison, the output from [JASP](https://jasp-stats.org/) on the 
    !! same data set is as follows.
    !! @image html Single_Factor_ANOVA_JASP_Output.png
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Analysis_of_variance)
    !! - [SPC Excel Single Factor ANOVA](https://www.spcforexcel.com/knowledge/root-cause-analysis/single-factor-anova)
    !! - [SPC Excel Gage R&R](https://www.spcforexcel.com/knowledge/measurement-systems-analysis/anova-gage-rr-part-1)
    !! - [SPC Excel Understanding Regression Statistics](https://www.spcforexcel.com/knowledge/root-cause-analysis/understanding-regression-statistics-part-1)
    !! - [NIST - Two Way ANOVA](https://www.itl.nist.gov/div898/handbook/prc/section4/prc427.htm)
    interface anova
        module procedure :: anova_1_factor
        module procedure :: anova_2_factor
        module procedure :: anova_model_fit
    end interface

    !> Computes the confidence interval for the specified distribution.
    !!
    !! @par Syntax Option 1
    !! @code{.f90}
    !! real(real64) function confidence_interval(class(distribution) dist, real(real64) alpha, real(real64) s, integer(int32) n)
    !! real(real32) function confidence_interval(class(distribution) dist, real(real32) alpha, real(real64) s, integer(int32) n)
    !! @endcode
    !!
    !! @param[in] dist The @ref distribution object defining the probability
    !! distribution to establish the confidence level.
    !! @param[in] alpha The probability value of interest.
    !! @param[in] s The sample standard deviation.
    !! @param[in] n The number of samples in the data set.
    !! @return The confidence interval.
    !!
    !! @par Syntax Option 2
    !! @par Syntax
    !! @code{.f90}
    !! real(real64) function confidence_interval(class(distribution) dist, real(real64) alpha, real(real64) x(:))
    !! real(real32) function confidence_interval(class(distribution) dist, real(real32) alpha, real(real32) x(:))
    !! @endcode
    !!
    !! @param[in] dist The @ref distribution object defining the probability
    !! distribution to establish the confidence level.
    !! @param[in] alpha The probability value of interest.
    !! @param[in] x An N-element array containing the data to analyze.
    !! @return The confidence interval.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Confidence_interval)
    interface confidence_interval
        module procedure :: confidence_interval_real64
        module procedure :: confidence_interval_real32
        module procedure :: confidence_interval_real64_array
        module procedure :: confidence_interval_real32_array
    end interface

    ! statistics_implementation.f90
    interface
        pure module function mean_real64(x) result(rst)
            real(real64), intent(in) :: x(:)
            real(real64) :: rst
        end function

        pure module function mean_real32(x) result(rst)
            real(real32), intent(in) :: x(:)
            real(real32) :: rst
        end function

        pure module function variance_real64(x) result(rst)
            real(real64), intent(in) :: x(:)
            real(real64) :: rst
        end function

        pure module function variance_real32(x) result(rst)
            real(real32), intent(in) :: x(:)
            real(real32) :: rst
        end function

        pure module function standard_deviation_real64(x) result(rst)
            real(real64), intent(in) :: x(:)
            real(real64) :: rst
        end function

        pure module function standard_deviation_real32(x) result(rst)
            real(real32), intent(in) :: x(:)
            real(real32) :: rst
        end function

        module function median_real64(x) result(rst)
            real(real64), intent(inout) :: x(:)
            real(real64) :: rst
        end function

        module function median_real32(x) result(rst)
            real(real32), intent(inout) :: x(:)
            real(real32) :: rst
        end function

        module function r_squared_real64(x, xm, err) result(rst)
            real(real64), intent(in) :: x(:), xm(:)
            class(errors), intent(inout), optional, target :: err
            real(real64) :: rst
        end function

        module function r_squared_real32(x, xm, err) result(rst)
            real(real32), intent(in) :: x(:), xm(:)
            class(errors), intent(inout), optional, target :: err
            real(real32) :: rst
        end function

        module function adjusted_r_squared_real64(p, x, xm, err) result(rst)
            integer(int32), intent(in) :: p
            real(real64), intent(in) :: x(:), xm(:)
            class(errors), intent(inout), optional, target :: err
            real(real64) :: rst
        end function

        module function adjusted_r_squared_real32(p, x, xm, err) result(rst)
            integer(int32), intent(in) :: p
            real(real32), intent(in) :: x(:), xm(:)
            class(errors), intent(inout), optional, target :: err
        end function

        pure module function quantile_real64(x, q) result(rst)
            real(real64), intent(in) :: x(:), q
            real(real64) :: rst
        end function

        pure module function quantile_real32(x, q) result(rst)
            real(real32), intent(in) :: x(:), q
            real(real32) :: rst
        end function

        pure module function confidence_interval_real64(dist, alpha, s, n) &
            result(rst)
            class(distribution), intent(in) :: dist
            real(real64), intent(in) :: alpha, s
            integer(int32), intent(in) :: n
            real(real64) :: rst
        end function

        pure module function confidence_interval_real32(dist, alpha, s, n) &
            result(rst)
            class(distribution), intent(in) :: dist
            real(real32), intent(in) :: alpha, s
            integer(int32), intent(in) :: n
            real(real32) :: rst
        end function

        pure module function confidence_interval_real64_array(dist, alpha, x) &
            result(rst)
            class(distribution), intent(in) :: dist
            real(real64), intent(in) :: alpha, x(:)
            real(real64) :: rst
        end function

        pure module function confidence_interval_real32_array(dist, alpha, x) &
            result(rst)
            class(distribution), intent(in) :: dist
            real(real32), intent(in) :: alpha, x(:)
            real(real32) :: rst
        end function
        
        module subroutine t_test_equal_var_real64(x1, x2, stat, p, dof)
            real(real64), intent(in) :: x1(:), x2(:)
            real(real64), intent(out) :: stat, p, dof
        end subroutine

        module subroutine t_test_equal_var_real32(x1, x2, stat, p, dof)
            real(real32), intent(in) :: x1(:), x2(:)
            real(real32), intent(out) :: stat, p, dof
        end subroutine

        module subroutine t_test_unequal_var_real64(x1, x2, stat, p, dof)
            real(real64), intent(in) :: x1(:), x2(:)
            real(real64), intent(out) :: stat, p, dof
        end subroutine

        module subroutine t_test_unequal_var_real32(x1, x2, stat, p, dof)
            real(real32), intent(in) :: x1(:), x2(:)
            real(real32), intent(out) :: stat, p, dof
        end subroutine

        module subroutine t_test_paired_real64(x1, x2, stat, p, dof, err)
            real(real64), intent(in) :: x1(:), x2(:)
            real(real64), intent(out) :: stat, p, dof
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module subroutine t_test_paired_real32(x1, x2, stat, p, dof, err)
            real(real32), intent(in) :: x1(:), x2(:)
            real(real32), intent(out) :: stat, p, dof
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module subroutine f_test_real64(x1, x2, stat, p, dof1, dof2)
            real(real64), intent(in) :: x1(:), x2(:)
            real(real64), intent(out) :: stat, p, dof1, dof2
        end subroutine

        module subroutine f_test_real32(x1, x2, stat, p, dof1, dof2)
            real(real32), intent(in) :: x1(:), x2(:)
            real(real32), intent(out) :: stat, p, dof1, dof2
        end subroutine
        
        module function anova_1_factor(x) result(rst)
            real(real64), intent(in) :: x(:,:)
            type(single_factor_anova_table) :: rst
        end function

        module function anova_2_factor(x) result(rst)
            real(real64), intent(in) :: x(:,:,:)
            type(two_factor_anova_table) :: rst
        end function

        module function anova_model_fit(nmodelparams, ymeas, ymod, err) result(rst)
            integer(int32), intent(in) :: nmodelparams
            real(real64), intent(in) :: ymeas(:), ymod(:)
            class(errors), intent(inout), optional, target :: err
            type(single_factor_anova_table) :: rst
        end function
    end interface

! ******************************************************************************
! SPECIAL FUNCTIONS
! ------------------------------------------------------------------------------
    !> Computes the beta function.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! real(real64) function beta(real(real64) a, real(real64) b)
    !! real(real32) function beta(real(real32) a, real(real32) b)
    !! @endcode
    !!
    !! @param[in] a The first argument of the function.
    !! @param[in] b The second argument of the function.
    !!
    !! @return The value of the beta function at @p a and @p b.
    !!
    !! @remarks The beta function is related to the gamma function
    !! by the following relationship \f$ \beta(a,b) = 
    !! \frac{\Gamma(a) \Gamma(b)}{\Gamma(a + b)} \f$.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Beta_function)
    interface beta
        module procedure :: beta_real64
        module procedure :: beta_real32
    end interface

    !> Computes the regularized beta function.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! real(real64) function regularized_beta(real(real64) a, real(real64) b, real(real64) x)
    !! real(real32) function regularized_beta(real(real32) a, real(real32) b, real(real32) x)
    !! @endcode
    !!
    !! @param[in] a The first argument of the function.
    !! @param[in] b The second argument of the function.
    !! @param[in] x The upper limit of the integration.
    !!
    !! @return The value of the incomplete beta function.
    !!
    !! @remarks The regularized beta function is defined as the ratio between
    !! the incomplete beta function and the beta function: \f$ I_{x}(a,b) = 
    !! \frac{\beta(x;a,b)}{\beta(a,b)} \f$.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Beta_function)
    interface regularized_beta
        module procedure :: regularized_beta_real64
        module procedure :: regularized_beta_real32
    end interface

    !> Computes the incomplete beta function.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! real(real64) function incomplete_beta(real(real64) a, real(real64) b, real(real64) x)
    !! real(real32) function incomplete_beta(real(real32) a, real(real32) b, real(real32) x)
    !! @endcode
    !!
    !! @param[in] a The first argument of the function.
    !! @param[in] b The second argument of the function.
    !! @param[in] x The upper limit of the integration.
    !!
    !! @return The value of the incomplete beta function.
    !!
    !! @remarks The incomplete beta function is defind as \f$ \beta(x;a,b) =
    !! \int_{0}^{x} t^{a-1} (1 - t)^{b-1} dt \f$.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Beta_function#Incomplete_beta_function)
    interface incomplete_beta
        module procedure :: incomplete_beta_real64
        module procedure :: incomplete_beta_real32
    end interface

    !> Computes the digamma function \f$ \psi(x) = 
    !! \frac{d}{dx}\left( \ln \left( \Gamma \left( x \right) \right) 
    !! \right) \f$.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! real(real64) function digamma(real(real64) x)
    !! real(real32) function digamma(real(real32) x)
    !! @endcode
    !!
    !! @param[in] x The value at which to evaluate the function.
    !! @return The function value.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Digamma_function)
    interface digamma
        module procedure :: digamma_real64
        module procedure :: digamma_real32
    end interface
    
    !> Computes the "upper" incomplete gamma function 
    !! \f$ \Gamma(a, x) = \int_{x}^{\infty} t^{a-1} e^{-t} \,dt \f$.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! real(real64) incomplete_gamma_upper(real(real64) a, real(real64) x)
    !! real(real32) incomplete_gamma_upper(real(real32) a, real(real32) x)
    !! @endcode
    !!
    !! @param[in] a The coefficient value.
    !! @param[in] x The value at which to evaluate the function.
    !! @return The function value.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Incomplete_gamma_function)
    interface incomplete_gamma_upper
        module procedure :: incomplete_gamma_upper_real64
        module procedure :: incomplete_gamma_upper_real32
    end interface

    !> Computes the "lower" incomplete gamma function 
    !! \f$ \gamma(a, x) = \int_{0}^{x} t^{a-1} e^{-t} \,dt \f$.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! real(real64) incomplete_gamma_lower(real(real64) a, real(real64) x)
    !! real(real32) incomplete_gamma_lower(real(real32) a, real(real32) x)
    !! @endcode
    !!
    !! @param[in] a The coefficient value.
    !! @param[in] x The value at which to evaluate the function.
    !! @return The function value.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Incomplete_gamma_function)
    interface incomplete_gamma_lower
        module procedure :: incomplete_gamma_lower_real64
        module procedure :: incomplete_gamma_lower_real32
    end interface

    ! special_functions_beta.f90
    interface
        pure elemental module function beta_real64(a, b) result(rst)
            real(real64), intent(in) :: a, b
            real(real64) :: rst
        end function

        pure elemental module function beta_real32(a, b) result(rst)
            real(real32), intent(in) :: a, b
            real(real32) :: rst
        end function

        pure elemental module function regularized_beta_real64(a, b, x) result(rst)
            real(real64), intent(in) :: a, b, x
            real(real64) :: rst
        end function

        pure elemental module function regularized_beta_real32(a, b, x) result(rst)
            real(real32), intent(in) :: a, b, x
            real(real32) :: rst
        end function

        pure elemental module function incomplete_beta_real64(a, b, x) result(rst)
            real(real64), intent(in) :: a, b, x
            real(real64) :: rst
        end function

        pure elemental module function incomplete_beta_real32(a, b, x) result(rst)
            real(real32), intent(in) :: a, b, x
            real(real32) :: rst
        end function
    end interface

    ! special_functions_digamma.f90
    interface
        pure elemental module function digamma_real64(x) result(rst)
            real(real64), intent(in) :: x
            real(real64) :: rst
        end function

        pure elemental module function digamma_real32(x) result(rst)
            real(real32), intent(in) :: x
            real(real32) :: rst
        end function
    end interface

    ! special_functions_gamma.f90
    interface
        pure elemental module function incomplete_gamma_upper_real64(a, x) &
            result(rst)
            real(real64), intent(in) :: a, x
            real(real64) :: rst
        end function

        pure elemental module function incomplete_gamma_upper_real32(a, x) &
            result(rst)
            real(real32), intent(in) :: a, x
            real(real32) :: rst
        end function

        pure elemental module function incomplete_gamma_lower_real64(a, x) &
            result(rst)
            real(real64), intent(in) :: a, x
            real(real64) :: rst
        end function

        pure elemental module function incomplete_gamma_lower_real32(a, x) &
            result(rst)
            real(real32), intent(in) :: a, x
            real(real32) :: rst
        end function
    end interface

! ******************************************************************************
! REGRESSION
! ------------------------------------------------------------------------------
    !> A container for regression-related statistical information.
    type regression_statistics
        !> The standard error for the model coefficient.
        !!
        !! \f$ se(\beta_{i}) = \sqrt{\sigma^{2} C_{ii}} \f$
        real(real64) :: standard_error
        !> The T-statistic for the model coefficient.
        !!
        !! \f$ t_o = \frac{ \beta_{i} }{se(\beta_{i})} \f$
        real(real64) :: t_statistic
        !> The probability that the coefficient is not statistically important.
        !! A statistically important coefficient will have a low probability 
        !! (p-value), typically 0.05 or lower; however, a p-value of up to ~0.2
        !! may be acceptable dependent upon the problem.  Typically any p-value
        !! larger than ~0.2 indicates the parameter is not statistically
        !! important for the model.
        !!
        !! \f$ p = t_{|t_o|, df_{residual}} \f$
        real(real64) :: probability
        !> The confidence interval for the parameter at the level determined by
        !! the regression process.
        !!
        !! \f$ c = t_{\alpha, df} se(\beta_{i}) \f$
        real(real64) :: confidence_interval
    end type

    !> Computes the coefficient matrix \f$ X \f$ to the linear least-squares 
    !! regression problem of \f$ X \beta = y \f$, where \f$ X \f$ is the 
    !! coefficient matrix computed here, \f$ \beta \f$ is the vector of 
    !! coefficients to be determined, and \f$ y \f$ is the vector of measured 
    !! dependent variables.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! subroutine coefficient_matrix(integer(int32) order, logical intercept, real(real64) x(:), real(real64) c(:,:), optional class(errors) err)
    !! subroutine coefficient_matrix(integer(int32) order, logical intercept, real(real32) x(:), real(real32) c(:,:), optional class(errors) err)
    !! @endcode
    !!
    !! @param[in] order The order of the equation to fit.  This value must be
    !!  at least one (linear equation), but can be higher as desired.
    !! @param[in] intercept Set to true if the intercept is being computed
    !!  as part of the regression; else, false.
    !! @param[in] x An N-element array containing the independent variable
    !!  measurement points.
    !! @param[out] c An N-by-K matrix where the results will be written.  K
    !!  must equal order + 1 in the event @p intercept is true; however, if
    !!  @p intercept is false, K must equal order.
    !! @param[in,out] err A mechanism for communicating errors and warnings
    !!  to the caller.  Possible warning and error codes are as follows.
    !! - FS_NO_ERROR: No errors encountered.
    !! - FS_ARRAY_SIZE_ERROR: Occurs if @p c is not properly sized.
    !! - FS_INVALID_INPUT_ERROR: Occurs if @p order is less than 1.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Linear_regression)
    interface coefficient_matrix
        module procedure :: coefficient_matrix_real64
        module procedure :: coefficient_matrix_real32
    end interface

    !> Computes the covariance matrix C where 
    !! \f$ C = \left( X^{T} X \right)^{-1} \f$ where \f$ X \f$ is computed
    !! by @ref coefficient_matrix.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! subroutine covariance_matrix(real(real64) x(:,:), real(real64) c(:,:), optional class(errors) err)
    !! subroutine covariance_matrix(real(real32) x(:,:), real(real32) c(:,:), optional class(errors) err)
    !! @endcode
    !!
    !! @param[in] x An M-by-N matrix containing the formatted independent data
    !!  matrix \f$ X \f$ as computed by @ref coefficient_matrix.
    !! @param[out] c The N-by-N covariance matrix.
    !! @param[in,out] err A mechanism for communicating errors and warnings
    !!  to the caller.  Possible warning and error codes are as follows.
    !! - FS_NO_ERROR: No errors encountered.
    !! - FS_ARRAY_SIZE_ERROR: Occurs if any of the matrices are not sized
    !!      correctly.
    !! - FS_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Covariance_matrix)
    !! - [Wikipedia - Regression](https://en.wikipedia.org/wiki/Linear_regression)
    interface covariance_matrix
        module procedure :: covariance_matrix_real64
        module procedure :: covariance_matrix_real32
    end interface

    !> Computes a linear least-squares regression to fit a set of data.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! subroutine linear_least_squares(integer(int32) order, lobical intercept, real(real64) x(:), real(real64) y(:), real(real64) coeffs(:), real(real64) ymod(:), real(real64) resid, optional type(regression_statistics) stats(:), optional real(real64) alpha, optional class(errors) err)
    !! subroutine linear_least_squares(integer(int32) order, lobical intercept, real(real32) x(:), real(real32) y(:), real(real32) coeffs(:), real(real32) ymod(:), real(real32) resid, optional type(regression_statistics) stats(:), optional real(real32) alpha, optional class(errors) err)
    !! @endcode
    !!
    !! @param[in] order The order of the equation to fit.  This value must be
    !!  at least one (linear equation), but can be higher as desired.
    !! @param[in] intercept Set to true if the intercept is being computed
    !!  as part of the regression; else, false.
    !! @param[in] x An N-element array containing the independent variable
    !!  measurement points.
    !! @param[in] y An N-element array containing the dependent variable
    !!  measurement points.
    !! @param[out] coeffs An ORDER+1 element array where the coefficients will
    !!  be written.
    !! @param[out] ymod An N-element array where the modeled data will be 
    !!  written.
    !! @param[out] resid An N-element array where the residual error data
    !!  will be written (modeled - actual).
    !! @param[out] stats An M-element array of @ref regression_statistics items
    !!  where M = ORDER + 1 when @p intercept is set to true; however, if 
    !!  @p intercept is set to false, M = ORDER.
    !! @param[in] alpha The significance level at which to evaluate the
    !!  confidence intervals.  The default value is 0.05 such that a 95% 
    !!  confidence interval is calculated.
    !! @param[in,out] err A mechanism for communicating errors and warnings
    !!  to the caller.  Possible warning and error codes are as follows.
    !! - FS_NO_ERROR: No errors encountered.
    !! - FS_ARRAY_SIZE_ERROR: Occurs if any of the arrays are not approriately
    !!      sized.
    !! - FS_INVALID_INPUT_ERROR: Occurs if @p order is less than 1.
    !! - FS_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Linear_regression)
    !! - [SPC Excel Understanding Regression Statistics](https://www.spcforexcel.com/knowledge/root-cause-analysis/understanding-regression-statistics-part-1)
    interface linear_least_squares
        module procedure :: linear_least_squares_real64
        module procedure :: linear_least_squares_real32
    end interface

    ! regression_implementation.f90
    interface
        module subroutine coefficient_matrix_real64(order, intercept, x, c, err)
            integer(int32), intent(in) :: order
            logical, intent(in) :: intercept
            real(real64), intent(in) :: x(:)
            real(real64), intent(out) :: c(:,:)
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module subroutine coefficient_matrix_real32(order, intercept, x, c, err)
            integer(int32), intent(in) :: order
            logical, intent(in) :: intercept
            real(real32), intent(in) :: x(:)
            real(real32), intent(out) :: c(:,:)
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module subroutine covariance_matrix_real64(x, c, err)
            real(real64), intent(in) :: x(:,:)
            real(real64), intent(out) :: c(:,:)
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module subroutine covariance_matrix_real32(x, c, err)
            real(real32), intent(in) :: x(:,:)
            real(real32), intent(out) :: c(:,:)
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module subroutine linear_least_squares_real64(order, intercept, x, y, &
            coeffs, ymod, resid, stats, alpha, err)
            integer(int32), intent(in) :: order
            logical, intent(in) :: intercept
            real(real64), intent(in) :: x(:), y(:)
            real(real64), intent(out) :: coeffs(:), ymod(:), resid(:)
            type(regression_statistics), intent(out), optional :: stats(:)
            real(real64), intent(in), optional :: alpha
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module subroutine linear_least_squares_real32(order, intercept, x, y, &
            coeffs, ymod, resid, stats, alpha, err)
            integer(int32), intent(in) :: order
            logical, intent(in) :: intercept
            real(real32), intent(in) :: x(:), y(:)
            real(real32), intent(out) :: coeffs(:), ymod(:), resid(:)
            type(regression_statistics), intent(out), optional :: stats(:)
            real(real32), intent(in), optional :: alpha
            class(errors), intent(inout), optional, target :: err
        end subroutine
    end interface

! ******************************************************************************
! EXPERIMENTAL DESIGN
! ------------------------------------------------------------------------------
    !> Computes the appropriate size for a full-factorial design table.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! subroutine get_full_factorial_matrix_size(integer(int32) vars(:), integer(int32) m, integer(int32) n, optional class(errors) err)
    !! @endcode
    !!
    !! @param[in] vars An M-element array containing the M factors to study.  
    !! Each of the M entries to the array is expected to contain the number of
    !! options for that particular factor to explore.  This value must be 
    !! greater than or equal to 1.
    !! @param[out] m The number of rows for the table.
    !! @param[out] n The number of columns for the table.
    !! @param[in,out] err A mechanism for communicating errors and warnings
    !!  to the caller.  Possible warning and error codes are as follows.
    !! - FS_NO_ERROR: No errors encountered.
    !! - FS_INVALID_INPUT_ERROR: Occurs if any items in @p vars are less than 1.
    interface get_full_factorial_matrix_size
        module procedure :: get_full_factorial_matrix_size_int32
    end interface

    !> Computes a table with values scaled from 1 to N describing a 
    !! full-factorial design.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! subroutine full_factorial(integer(int32) vars(:), integer(int32) tbl(:,:), optional class(errors) err)
    !! @endcode
    !!
    !! @param[in] vars An M-element array containing the M factors to study.  
    !! Each of the M entries to the array is expected to contain the number of
    !! options for that particular factor to explore.  This value must be 
    !! greater than or equal to 1.
    !! @param[out] tbl A table where the design will be written.  Use 
    !! @ref get_full_factorial_matrix_size to determine the appropriate table
    !! size.
    !! @param[in,out] err A mechanism for communicating errors and warnings
    !!  to the caller.  Possible warning and error codes are as follows.
    !! - FS_NO_ERROR: No errors encountered.
    !! - FS_INVALID_INPUT_ERROR: Occurs if any items in @p vars are less than 1.
    !! - FS_ARRAY_SIZE_ERROR: Occurs if @p tbl is not properly sized.
    !!
    !! @par Example
    !! The following example illustrates how to construct a full-factorial 
    !! design consisting of 3 distinct factors with different design options for
    !! each factor.
    !! @code{.f90}
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
    !! @endcode
    !! The above program produces the following output.
    !! @code{.txt}
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
    !! @endcode
    !! The resulting table is simply integer values representing each design
    !! point.  It is up to the user to associate these integer values with
    !! meaningful design parameters.
    interface full_factorial
        module procedure :: full_factorial_int32
    end interface

    ! experimental_design_implementation.f90
    interface
        module subroutine get_full_factorial_matrix_size_int32(vars, m, n, err)
            integer(int32), intent(in) :: vars(:)
            integer(int32), intent(out) :: m, n
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module subroutine full_factorial_int32(vars, tbl, err)
            integer(int32), intent(in) :: vars(:)
            integer(int32), intent(out) :: tbl(:,:)
            class(errors), intent(inout), optional, target :: err
        end subroutine
    end interface

! ------------------------------------------------------------------------------
end module
! fstats.f90

!> @mainpage
!!
!! @section intro_sec Introduction
!! FSTATS is a modern Fortran statistical library containing routines for 
!! computing basic statistical properties, hypothesis testing, regression, 
!! special functions, and even experimental design.
!!
!! @par Regression Example
!! The following example illustrates fitting a cubic polynomial to a set 
!! of data points.
!! @code{.f90}
!! program example
!!     use iso_fortran_env
!!     use fstats
!!     implicit none
!!
!!     ! Local Variables
!!     character, parameter :: tab = achar(9)
!!     character, parameter :: nl = new_line('a')
!!     integer(int32) :: i
!!     real(real64) :: x(31), y(31), coeffs(4), ymodeled(31), residuals(31)
!!     type(regression_statistics) :: stats(4)
!!
!!     ! Define the data
!!     x = [ &
!!             0.0d0, &
!!             0.1d0, &
!!             0.2d0, &
!!             0.3d0, &
!!             0.4d0, &
!!             0.5d0, &
!!             0.6d0, &
!!             0.7d0, &
!!             0.8d0, &
!!             0.9d0, &
!!             1.0d0, &
!!             1.1d0, &
!!             1.2d0, &
!!             1.3d0, &
!!             1.4d0, &
!!             1.5d0, &
!!             1.6d0, &
!!             1.7d0, &
!!             1.8d0, &
!!             1.9d0, &
!!             2.0d0, &
!!             2.1d0, &
!!             2.2d0, &
!!             2.3d0, &
!!             2.4d0, &
!!             2.5d0, &
!!             2.6d0, &
!!             2.7d0, &
!!             2.8d0, &
!!             2.9d0, &
!!             3.0d0 &
!!         ]
!!         y = [ &
!!             0.577855138260449d0, &
!!             0.614883095604222d0, &
!!             0.633891127488559d0, &
!!             0.718405829701721d0, &
!!             0.753668502759107d0, &
!!             0.814967857310145d0, &
!!             0.861870996499704d0, &
!!             0.925100533744381d0, &
!!             0.947038018520063d0, &
!!             1.025198043343280d0, &
!!             1.042142354497610d0, &
!!             1.121528566784440d0, &
!!             1.177570314994070d0, &
!!             1.229237567525370d0, &
!!             1.261114062593870d0, &
!!             1.296408162551430d0, &
!!             1.394353657051120d0, &
!!             1.367144391560370d0, &
!!             1.428164431435150d0, &
!!             1.548944935073270d0, &
!!             1.505100149282990d0, &
!!             1.560701023751520d0, &
!!             1.609113012481530d0, &
!!             1.663687366875500d0, &
!!             1.707149545456870d0, &
!!             1.800935947618110d0, &
!!             1.847819988906440d0, &
!!             1.884242821675810d0, &
!!             1.966174239373140d0, &
!!             1.977005266443110d0, &
!!             2.034137257154140d0 &    
!!         ]
!!
!!         ! Fit the data
!!         call linear_least_squares(3, .true., x, y, coeffs, ymodeled, residuals, stats)
!!
!!         ! Display the results
!!         print '(AF8.5AF8.5AF8.5AF8.5A)', "Model: y = ", &
!!             coeffs(1), " + ", &
!!             coeffs(2), " x + ", &
!!             coeffs(3), " x**2 + ", &
!!             coeffs(4), " x**3"
!!
!!         ! Illustrate the statistics for each coefficient
!!         do i = 1, size(stats)
!!             print '(AI0AF6.3AF6.3AF6.3AF6.3)', &
!!                 "Coefficient ", i, ":" // nl // &
!!                 tab // "Standard Error: ", stats(i)%standard_error, nl // &
!!                 tab // "Confidence Interval: +/-", stats(i)%confidence_interval, nl // &
!!                 tab // "T-Statistic: ", stats(i)%t_statistic, nl // &
!!                 tab // "P-Value: ", stats(i)%probability
!!         end do
!! end program
!! @endcode
!! The above program produces the following results.
!! @code{.txt}
!! Model: y =  0.55341 +  0.55399 x + -0.05483 x**2 +  0.01185 x**3
!! Coefficient 1:
!!         Standard Error:  0.015
!!         Confidence Interval: +/- 0.031
!!         T-Statistic: 36.404
!!         P-Value:  0.000
!! Coefficient 2:
!!         Standard Error:  0.045
!!         Confidence Interval: +/- 0.092
!!         T-Statistic: 12.417
!!         P-Value:  0.000
!! Coefficient 3:
!!         Standard Error:  0.035
!!         Confidence Interval: +/- 0.072
!!         T-Statistic: -1.572
!!         P-Value:  0.128
!! Coefficient 4:
!!         Standard Error:  0.008
!!         Confidence Interval: +/- 0.016
!!         T-Statistic:  1.552
!!         P-Value:  0.132
!! @endcode
!! Plotting the data results in the following plot.
!! @image html Polynomial_Fit_Results.png"

!> @brief Provides various statistical routines.
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
    public :: chi_squared_distribution
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
    public :: FS_NO_ERROR
    public :: FS_ARRAY_SIZE_ERROR
    public :: FS_INVALID_INPUT_ERROR
    public :: FS_OUT_OF_MEMORY_ERROR
    public :: FS_UNDERDEFINED_PROBLEM_ERROR
    public :: FS_TOLERANCE_TOO_SMALL_ERROR
    public :: FS_TOO_FEW_ITERATION_ERROR
    public :: FS_LEVENBERG_MARQUARDT_UPDATE
    public :: FS_QUADRATIC_UPDATE
    public :: FS_NIELSEN_UPDATE
    public :: assignment(=)

! ******************************************************************************
! ERROR CODES
! ------------------------------------------------------------------------------
    integer(int32), parameter :: FS_NO_ERROR = 0
    integer(int32), parameter :: FS_ARRAY_SIZE_ERROR = 10000
    integer(int32), parameter :: FS_INVALID_INPUT_ERROR = 10001
    integer(int32), parameter :: FS_OUT_OF_MEMORY_ERROR = 10002
    integer(int32), parameter :: FS_UNDERDEFINED_PROBLEM_ERROR = 10003
    integer(int32), parameter :: FS_TOLERANCE_TOO_SMALL_ERROR = 10004
    integer(int32), parameter :: FS_TOO_FEW_ITERATION_ERROR = 10005

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

    !> @brief Defines a Chi-squared distribution.
    !!
    !! @par Remarks
    !! For additional information see the following references.
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Chi-squared_distribution)
    type, extends(distribution) :: chi_squared_distribution
        !> The number of degrees of freedom.
        integer(int32) :: dof
    contains
        !> Computes the probability density function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function pdf(class(chi_squared_distribution) this, real(real64 )x)
        !! @endcode
        !! 
        !! @param[in] this The chi_squared_distribution object.
        !! @param[in] x The value at which to evaluate the function.
        !! @return The value of the function.
        !!
        !! @remarks The PDF for a Chi-squared distribution is given as 
        !! \f$ f(x) = \frac{x^{k/2 - 1} \exp{-x / 2}} {2^{k / 2} 
        !! \Gamma \left( \frac{k}{2} \right)} \f$, where \f$ k \f$ is the 
        !! number of degrees of freedom.
        procedure, public :: pdf => cs_pdf
        !> Computes the cumulative distribution function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function cdf(class(chi_squared_distribution) this, real(real64 )x)
        !! @endcode
        !! 
        !! @param[in] this The chi_squared_distribution object.
        !! @param[in] x The value at which to evaluate the function.
        !! @return The value of the function.
        !!
        !! @remarks The CDF for a Chi-squared distribution is given as 
        !! \f$ F(x) = \frac{ \gamma \left( \frac{k}{2}, \frac{x}{2} \right) }
        !! { \Gamma \left( \frac{k}{2} \right)} \f$
        !! where \f$ \gamma \f$ is the 
        !! <a href = "https://en.wikipedia.org/wiki/Incomplete_gamma_function">
        !! lower incomplete gamma function</a>, and \f$ k \f$ is the number of 
        !! degrees of freedom.
        procedure, public :: cdf => cs_cdf
        !> Computes the mean of the distribution.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function mean(class(chi_squared_distribution) this)
        !! @endcode
        !!
        !! @param[in] this The chi_squared_distribution object.
        !! @return The mean value.
        procedure, public :: mean => cs_mean
        !> Computes the median of the distribution.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function median(class(chi_squared_distribution) this)
        !! @endcode
        !!
        !! @param[in] this The chi_squared_distribution object.
        !! @return The median value.
        procedure, public :: median => cs_median
        !> Computes the mode of the distribution.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function mode(class(chi_squared_distribution) this)
        !! @endcode
        !!
        !! @param[in] this The chi_squared_distribution object.
        !! @return The mode value.
        procedure, public :: mode => cs_mode
        !> Computes the variance of the distribution.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function variance(class(chi_squared_distribution) this)
        !! @endcode
        !!
        !! @param[in] this The chi_squared_distribution object.
        !! @return The variance value.
        procedure, public :: variance => cs_variance
    end type

    ! distributions_chisquared.f90
    interface
        pure module elemental function cs_pdf(this, x) result(rst)
            class(chi_squared_distribution), intent(in) :: this
            real(real64), intent(in) :: x
            real(real64) :: rst
        end function

        pure module elemental function cs_cdf(this, x) result(rst)
            class(chi_squared_distribution), intent(in) :: this
            real(real64), intent(in) :: x
            real(real64) :: rst
        end function

        pure module function cs_mean(this) result(rst)
            class(chi_squared_distribution), intent(in) :: this
            real(real64) :: rst
        end function

        pure module function cs_median(this) result(rst)
            class(chi_squared_distribution), intent(in) :: this
            real(real64) :: rst
        end function

        pure module function cs_mode(this) result(rst)
            class(chi_squared_distribution), intent(in) :: this
            real(real64) :: rst
        end function

        pure module function cs_variance(this) result(rst)
            class(chi_squared_distribution), intent(in) :: this
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
    !! The following example illustrates a single-factor ANOVA on the following
    !! data set.
    !! @image html Single_Factor_ANOVA_Data_Set.png
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
    !! subroutine linear_least_squares(integer(int32) order, logical intercept, real(real64) x(:), real(real64) y(:), real(real64) coeffs(:), real(real64) ymod(:), real(real64) resid, optional type(regression_statistics) stats(:), optional real(real64) alpha, optional class(errors) err)
    !! subroutine linear_least_squares(integer(int32) order, logical intercept, real(real32) x(:), real(real32) y(:), real(real32) coeffs(:), real(real32) ymod(:), real(real32) resid, optional type(regression_statistics) stats(:), optional real(real32) alpha, optional class(errors) err)
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
    !! @par Example
    !! The following example illustrates fitting a cubic polynomial to a set 
    !! of data points.
    !! @image html Polynomial_Data_Set.png
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use fstats
    !!     implicit none
    !!
    !!     ! Local Variables
    !!     character, parameter :: tab = achar(9)
    !!     character, parameter :: nl = new_line('a')
    !!     integer(int32) :: i
    !!     real(real64) :: x(31), y(31), coeffs(4), ymodeled(31), residuals(31)
    !!     type(regression_statistics) :: stats(4)
    !!
    !!     ! Define the data
    !!     x = [ &
    !!             0.0d0, &
    !!             0.1d0, &
    !!             0.2d0, &
    !!             0.3d0, &
    !!             0.4d0, &
    !!             0.5d0, &
    !!             0.6d0, &
    !!             0.7d0, &
    !!             0.8d0, &
    !!             0.9d0, &
    !!             1.0d0, &
    !!             1.1d0, &
    !!             1.2d0, &
    !!             1.3d0, &
    !!             1.4d0, &
    !!             1.5d0, &
    !!             1.6d0, &
    !!             1.7d0, &
    !!             1.8d0, &
    !!             1.9d0, &
    !!             2.0d0, &
    !!             2.1d0, &
    !!             2.2d0, &
    !!             2.3d0, &
    !!             2.4d0, &
    !!             2.5d0, &
    !!             2.6d0, &
    !!             2.7d0, &
    !!             2.8d0, &
    !!             2.9d0, &
    !!             3.0d0 &
    !!         ]
    !!         y = [ &
    !!             0.577855138260449d0, &
    !!             0.614883095604222d0, &
    !!             0.633891127488559d0, &
    !!             0.718405829701721d0, &
    !!             0.753668502759107d0, &
    !!             0.814967857310145d0, &
    !!             0.861870996499704d0, &
    !!             0.925100533744381d0, &
    !!             0.947038018520063d0, &
    !!             1.025198043343280d0, &
    !!             1.042142354497610d0, &
    !!             1.121528566784440d0, &
    !!             1.177570314994070d0, &
    !!             1.229237567525370d0, &
    !!             1.261114062593870d0, &
    !!             1.296408162551430d0, &
    !!             1.394353657051120d0, &
    !!             1.367144391560370d0, &
    !!             1.428164431435150d0, &
    !!             1.548944935073270d0, &
    !!             1.505100149282990d0, &
    !!             1.560701023751520d0, &
    !!             1.609113012481530d0, &
    !!             1.663687366875500d0, &
    !!             1.707149545456870d0, &
    !!             1.800935947618110d0, &
    !!             1.847819988906440d0, &
    !!             1.884242821675810d0, &
    !!             1.966174239373140d0, &
    !!             1.977005266443110d0, &
    !!             2.034137257154140d0 &    
    !!         ]
    !!
    !!         ! Fit the data
    !!         call linear_least_squares(3, .true., x, y, coeffs, ymodeled, residuals, stats)
    !!
    !!         ! Display the results
    !!         print '(AF8.5AF8.5AF8.5AF8.5A)', "Model: y = ", &
    !!             coeffs(1), " + ", &
    !!             coeffs(2), " x + ", &
    !!             coeffs(3), " x**2 + ", &
    !!             coeffs(4), " x**3"
    !!
    !!         ! Illustrate the statistics for each coefficient
    !!         do i = 1, size(stats)
    !!             print '(AI0AF6.3AF6.3AF6.3AF6.3)', &
    !!                 "Coefficient ", i, ":" // nl // &
    !!                 tab // "Standard Error: ", stats(i)%standard_error, nl // &
    !!                 tab // "Confidence Interval: +/-", stats(i)%confidence_interval, nl // &
    !!                 tab // "T-Statistic: ", stats(i)%t_statistic, nl // &
    !!                 tab // "P-Value: ", stats(i)%probability
    !!         end do
    !! end program
    !! @endcode
    !! The above program produces the following results.
    !! @code{.txt}
    !! Model: y =  0.55341 +  0.55399 x + -0.05483 x**2 +  0.01185 x**3
    !! Coefficient 1:
    !!         Standard Error:  0.015
    !!         Confidence Interval: +/- 0.031
    !!         T-Statistic: 36.404
    !!         P-Value:  0.000
    !! Coefficient 2:
    !!         Standard Error:  0.045
    !!         Confidence Interval: +/- 0.092
    !!         T-Statistic: 12.417
    !!         P-Value:  0.000
    !! Coefficient 3:
    !!         Standard Error:  0.035
    !!         Confidence Interval: +/- 0.072
    !!         T-Statistic: -1.572
    !!         P-Value:  0.128
    !! Coefficient 4:
    !!         Standard Error:  0.008
    !!         Confidence Interval: +/- 0.016
    !!         T-Statistic:  1.552
    !!         P-Value:  0.132
    !! @endcode
    !! Plotting the data results in the following plot.
    !! @image html Polynomial_Fit_Results.png"
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Linear_regression)
    !! - [SPC Excel Understanding Regression Statistics](https://www.spcforexcel.com/knowledge/root-cause-analysis/understanding-regression-statistics-part-1)
    interface linear_least_squares
        module procedure :: linear_least_squares_real64
        module procedure :: linear_least_squares_real32
    end interface

    !> @brief Computes statistics for the quality of fit for a regression model.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! type(regression_statistics) function calculate_regression_statistics(real(real64) resid(:), real(real64) params(:), real(real64) c(:,:), optional real(real64) alpha, optional class(errors) err)
    !! type(regression_statistics) function calculate_regression_statistics(real(real32) resid(:), real(real32) params(:), real(real32) c(:,:), optional real(real32) alpha, optional class(errors) err)
    !! @endcode
    !!
    !! @param[in] resid An M-element array containing the model residual errors.
    !! @param[in] params An N-element array containing the model parameters.
    !! @param[in] c The N-by-N covariance matrix.
    !! @param[in] alpha An optional input describing the probability level of 
    !!  the confidence interval analysis.
    !! @param[in,out] err A mechanism for communicating errors and warnings
    !!  to the caller.  Possible warning and error codes are as follows.
    !! - FS_NO_ERROR: No errors encountered.
    !! - FS_ARRAY_SIZE_ERROR: Occurs if @p c is not sized correctly.
    !! - FS_OUT_OF_MEMORY_ERROR: Occurs if a memory error is encountered.
    !! @return A @ref regression_statistics object containing the analysis
    !!  results.
    interface calculate_regression_statistics
        module procedure :: calculate_regression_stats_r64
        module procedure :: calculate_regression_stats_r32
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

        module function calculate_regression_stats_r64(resid, params, c, &
            alpha, err) result(rst)
            real(real64), intent(in) :: resid(:), params(:), c(:,:)
            real(real64), intent(in), optional :: alpha
            class(errors), intent(inout), optional, target :: err
            type(regression_statistics), allocatable :: rst(:)
        end function

        module function calculate_regression_stats_r32(resid, params, c, &
            alpha, err) result(rst)
            real(real32), intent(in) :: resid(:), params(:), c(:,:)
            real(real32), intent(in), optional :: alpha
            class(errors), intent(inout), optional, target :: err
            type(regression_statistics), allocatable :: rst(:)
        end function
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

! ******************************************************************************
! NONLINEAR REGRESSION
! ------------------------------------------------------------------------------
    !> @brief Provides a collection of iteration control parameters.
    type iteration_controls
        !> Defines the maximum number of iterations allowed.
        integer(int32) :: max_iteration_count
        !> Defines the maximum number of function evaluations allowed.
        integer(int32) :: max_function_evaluations
        !> Defines a tolerance on the gradient of the fitted function.
        real(real64) :: gradient_tolerance
        !> Defines a tolerance on the change in parameter values.
        real(real64) :: change_in_solution_tolerance
        !> Defines a tolerance on the metric associated with the residual error.
        real(real64) :: residual_tolerance
        !> Defines a tolerance to ensure adequate improvement on each iteration.
        real(real64) :: iteration_improvement_tolerance
        !> Defines how many iterations can pass before a re-evaluation of the
        !! Jacobian matrix is forced.
        integer(int32) :: max_iteration_between_updates
    contains
        !> @brief Sets the object to its default values.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_to_default(class(iteration_controls) x)
        !! @endcode
        !!
        !! @param[in,out] x The @ref iteration_controls object.
        procedure, public :: set_to_default => lm_set_default_tolerances
    end type

    !> @brief Provides information regarding convergence status.
    type convergence_info
        !> True if convergence on the gradient was achieved; else, false.
        logical :: converge_on_gradient
        !> The value of the gradient test parameter.
        real(real64) :: gradient_value
        !> True if convergence on the change in solution was achieved; else,
        !! false.
        logical :: converge_on_solution_change
        !> The value of the change in solution parameter.
        real(real64) :: solution_change_value
        !> True if convergence on the residual error parameter was achieved; 
        !! else, false.
        logical :: converge_on_residual_parameter
        !> The value of the residual error parameter.
        real(real64) :: residual_value
        !> True if the solution did not converge in the allowed number of 
        !! iterations.
        logical :: reach_iteration_limit
        !> The iteration count.
        integer(int32) :: iteration_count
        !> True if the solution did not converge in the allowed number of
        !! function evaluations.
        logical :: reach_function_evaluation_limit
        !> The function evaluation count.
        integer(int32) :: function_evaluation_count
        !> True if the user requested the stop; else, false.
        logical :: user_requested_stop
    end type

    !> @brief Options to control the Levenberg-Marquardt solver.
    type lm_solver_options
        !> The solver method to utilize.
        !! - FS_LEVENBERG_MARQUARDT_UPDATE:
        !! - FS_QUADRATIC_UPDATE:
        !! - FS_NIELSEN_UDPATE:
        integer(int32) :: method
        !> The step size used for the finite difference calculations of the
        !! Jacobian matrix.
        real(real64) :: finite_difference_step_size
        !> The factor to use when increasing the damping parameter.
        real(real64) :: damping_increase_factor
        !> The factor to use when decreasing the damping parameter.
        real(real64) :: damping_decrease_factor
    contains
        !> @brief Sets the object to its default values.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_to_default(class(lm_solver_options) x)
        !! @endcode
        !!
        !! @param[in,out] x The @ref lm_solver_options object.
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
            procedure(regression_function), intent(in), pointer :: fun
            real(real64), intent(in) :: xdata(:), params(:)
            real(real64), intent(out) :: jac(:,:)
            logical, intent(out) :: stop
            real(real64), intent(in), optional, target :: f0(:)
            real(real64), intent(out), optional, target :: f1(:)
            real(real64), intent(in), optional :: step
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module subroutine nonlinear_least_squares_1(fun, x, y, params, ymod, &
            resid, weights, maxp, minp, stats, alpha, controls, settings, &
            info, status, err)
            procedure(regression_function), intent(in), pointer :: fun
            real(real64), intent(in) :: x(:), y(:)
            real(real64), intent(inout) :: params(:)
            real(real64), intent(out) :: ymod(:), resid(:)
            real(real64), intent(in), optional, target :: weights(:), maxp(:), &
                minp(:)
            type(regression_statistics), intent(out), optional :: stats(:)
            real(real64), intent(in), optional :: alpha
            type(iteration_controls), intent(in), optional :: controls
            type(lm_solver_options), intent(in), optional :: settings
            type(convergence_info), intent(out), optional, target :: info
            procedure(iteration_update), intent(in), pointer, optional :: status
            class(errors), intent(inout), optional, target :: err
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
            class(iteration_controls), intent(inout) :: x
        end subroutine

        module subroutine lm_set_default_settings(x)
            class(lm_solver_options), intent(inout) :: x
        end subroutine
    end interface

    !> @brief Computes the Jacobian matrix for a nonlinear regression problem.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! subroutine jacobian(pointer procedure(regression_function) fun, real(real64) xdata(:), real(real64) params(:), real(real64) jac(:,:), logical stop, optional real(real64) f0(:), optional real(real64) f1(:), optional real(real64) step, optional class(errors) err)
    !! @endcode
    !!
    !! @param[in] fun A pointer to the @ref regression_function to evaluate.
    !! @param[in] xdata The M-element array containing x-coordinate data.
    !! @param[in] params The N-element array containing the model parameters.
    !! @param[out] jac The M-by-N matrix where the Jacobian will be written.
    !! @param[out] stop A value that the user can set in @p fun forcing the
    !!  evaluation process to stop prior to completion.
    !! @param[in] f0 An optional M-element array containing the model values
    !!  using the current parameters as defined in @p m.  This input can be
    !!  used to prevent the routine from performing a function evaluation at the
    !!  model parameter state defined in @p params.
    !! @param[out] f1 An optional M-element workspace array used for function
    !!  evaluations.
    !! @param[in] step The differentiation step size.  The default is the
    !!  square root of machine precision.
    !! @param[in,out] err A mechanism for communicating errors and warnings
    !!  to the caller.  Possible warning and error codes are as follows.
    !! - FS_NO_ERROR: No errors encountered.
    !! - FS_ARRAY_SIZE_ERROR: Occurs if any of the arrays are not properly
    !!      sized.
    !! - FS_OUT_OF_MEMORY_ERROR: Occurs if a memory error is encountered.
    interface jacobian
        module procedure :: regression_jacobian_1
    end interface

    !> @brief Performs a nonlinear regression to fit a model using a version
    !! of the Levenberg-Marquardt algorithm.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! subroutine nonlinear_least_squares(pointer procedure(regression_function) fun, real(real64) x(:), real(real64) y(:), real(real64) params(:), real(real64) ymod(:), real(real64) resid(:), optional real(real64) weights(:), optional real(real64) maxp(:), optional real(real64) minp(:), optional type(regression_statistics) stats(:), optional real(real64) alpha, optional type(iteration_controls) controls, optional type(lm_solver_options) settings, optional type(convergence_info) info, optional pointer procedure(iteration_update) status, optional class(errors) err)
    !! @endcode
    !!
    !! @param[in] fun A pointer to the @ref regression_function to evaluate.
    !! @param[in] x The M-element array containing x-coordinate data.
    !! @param[in] y The M-element array containing y-coordinate data.
    !! @param[in,out] params On input, the N-element array containing the
    !!  initial estimate of the model parameters.  On output, the computed 
    !!  model parameters.
    !! @param[out] ymod An M-element array where the modeled dependent data will
    !!  be written.
    !! @param[out] resid An M-element array where the model residuals will be
    !!  written.
    !! @param[in] weights An optional M-element array allowing the weighting of
    !!  individual points.
    !! @param[in] maxp An optional N-element array that can be used as upper
    !!  limits on the parameter values.  If no upper limit is requested for a
    !!  particular parameter, utilize a very large value.  The internal default
    !!  is to utilize huge() as a value.
    !! @param[in] minp An optional N-element array that can be used as lower
    !!  limits on the parameter values.  If no lower limit is requested for a
    !!  particalar parameter, utilize a very large magnitude, but negative,
    !!  value.  The internal default is to utilize -huge() as a value.
    !! @arapam[out] stats An optional N-element array that, if supplied, will
    !!  be used to return statistics about the fit for each parameter.
    !! @param[in] alpha An optional input describing the probability level of 
    !!  the confidence interval analysis assuming that statistics for each
    !!  parameter are being calculated.
    !! @param[in] controls An optional input providing custom iteration 
    !!  controls.
    !! @param[in] settings An optional input providing custom settings for the 
    !!  solver.
    !! @param[out] info An optional output that can be used to gain information
    !!  about the iterative solution and the nature of the convergence.
    !! @param[in] status An optional pointer to a routine that can be used to
    !!  extract iteration information.
    !! @param[in,out] err A mechanism for communicating errors and warnings
    !!  to the caller.  Possible warning and error codes are as follows.
    !! - FS_NO_ERROR: No errors encountered.
    !! - FS_ARRAY_SIZE_ERROR: Occurs if any of the arrays are not properly
    !!      sized.
    !! - FS_OUT_OF_MEMORY_ERROR: Occurs if a memory error is encountered.
    !! - FS_UNDERDEFINED_PROBLEM_ERROR: Occurs if the problem posed is 
    !!      underdetetermined (M < N).
    !! - FS_TOLERANCE_TOO_SMALL_ERROR: Occurs if any supplied tolerances are
    !!      too small to be practical.
    !! - FS_TOO_FEW_ITERATION_ERROR: Occurs if too few iterations are allowed.
    !!
    !! @par Example
    !! The following example illustrates the use of the nonlinear_least_squares
    !! routine for fitting a polynomial.  For comparison, the 
    !! @ref linear_least_squares routine is also employed.
    !!
    !! @code{.f90}
    !! module nl_example
    !!     use iso_fortran_env
    !! contains
    !!     subroutine exfun(x, p, f, stop)
    !!         ! Arguments
    !!         real(real64), intent(in) :: x(:), p(:)
    !!         real(real64), intent(out) :: f(:)
    !!         logical, intent(out) :: stop
    !!
    !!         ! Function
    !!         f = p(4) * x**3 + p(3) * x**2 + p(2) * x + p(1)
    !!
    !!         ! Do not stop
    !!         stop = .false.
    !!     end subroutine
    !! end module
    !!
    !! program example
    !!     use iso_fortran_env
    !!     use fstats
    !!     use nl_example
    !!     implicit none
    !!
    !!     ! Local Variables
    !!     character, parameter :: tab = achar(9)
    !!     character, parameter :: nl = new_line('a')
    !!     integer(int32) :: i
    !!     procedure(regression_function), pointer :: fun
    !!     real(real64) :: xp(21), yp(21), params(4), ymod(21), resid(21), &
    !!         ylmod(21), lresid(21), lparams(4)
    !!     type(convergence_info) :: info
    !!     type(regression_statistics) :: stats(4), lstats(4)
    !!
    !!     ! Data to fit
    !!     xp = [0.0d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.8d0, &
    !!         0.9d0, 1.0d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0, 1.5d0, 1.6d0, 1.7d0, &
    !!         1.8d0, 1.9d0, 2.0d0]
    !!     yp = [1.216737514d0, 1.250032542d0, 1.305579195d0, 1.040182335d0, &
    !!         1.751867738d0, 1.109716707d0, 2.018141531d0, 1.992418729d0, &
    !!         1.807916923d0, 2.078806005d0, 2.698801324d0, 2.644662712d0, &
    !!         3.412756702d0, 4.406137221d0, 4.567156645d0, 4.999550779d0, &
    !!         5.652854194d0, 6.784320119d0, 8.307936836d0, 8.395126494d0, &
    !!         10.30252404d0]
    !!
    !!     ! Assign the function pointer and define an initial solution estimate
    !!     fun => exfun
    !!     params = 1.0d0
    !!
    !!     ! Solve the problem
    !!     call nonlinear_least_squares(fun, xp, yp, params, ymod, resid, &
    !!         info = info, stats = stats)
    !!
    !!     ! Display the results
    !!     print '(A)', "Model:"
    !!     print 100, (tab // "params(", i, "): ", params(i), i = 1, size(params))
    !!
    !!     print '(A)', "Statistics:"
    !!     print 101, ( &
    !!         "Coefficient ", i, ":" // nl // &
    !!         tab // "Standard Error: ", stats(i)%standard_error, nl // &
    !!         tab // "Confidence Interval: +/-", stats(i)%confidence_interval, nl // &
    !!         tab // "T-Statistic: ", stats(i)%t_statistic, nl // &
    !!         tab // "P-Value: ", stats(i)%probability, &
    !!         i = 1, size(stats) &
    !!     )
    !!
    !!     print '(A)', "Convergence Information:"
    !!     print 102, tab // "Iteration Count: ", info%iteration_count
    !!     print 102, tab // "Function Evaluation Count: ", &
    !!         info%function_evaluation_count
    !!     print 103, tab // "Converge on Gradient: ", &
    !!         info%converge_on_gradient, " (", info%gradient_value, ")"
    !!     print 103, tab // "Converge on Solution Change: ", &
    !!         info%converge_on_solution_change, " (", info%solution_change_value, ")"
    !!     print 103, tab // "Converge on Parameter Change: ", &
    !!         info%converge_on_residual_parameter, " (", info%residual_value, ")"
    !!
    !!     ! As our model is simply a 3rd order polynomial, the linear_least_squares
    !!     ! routine can also be used.  As a comparison, here's the 
    !!     ! linear_least_squares solution
    !!     call linear_least_squares(3, .true., xp, yp, lparams, ylmod, lresid, lstats)
    !!
    !!     ! Print the linear model results
    !!     print '(A)', nl // "Linear Model:"
    !!     print 100, (tab // "params(", i, "): ", lparams(i), i = 1, size(lparams))
    !!
    !!     print '(A)', "Statistics:"
    !!     print 101, ( &
    !!         "Coefficient ", i, ":" // nl // &
    !!         tab // "Standard Error: ", lstats(i)%standard_error, nl // &
    !!         tab // "Confidence Interval: +/-", lstats(i)%confidence_interval, nl // &
    !!         tab // "T-Statistic: ", lstats(i)%t_statistic, nl // &
    !!         tab // "P-Value: ", lstats(i)%probability, &
    !!         i = 1, size(lstats) &
    !!     )
    !!
    !! 100 format(A, I0, A, F8.5)
    !! 101 format(A, I0, A, F6.3, A, F6.3, A, F6.3, A, F6.3)
    !! 102 format(A, I0)
    !! 103 format(A, L1, A, E9.3, A)
    !! end program
    !! @endcode
    !! The above program produces the following results.
    !! @code{.txt}
    !! Model:
    !!         params(1):  1.18661
    !!         params(2):  0.44661
    !!         params(3): -0.12232
    !!         params(4):  1.06476
    !! Statistics:
    !! Coefficient 1:
    !!         Standard Error:  0.230
    !!         Confidence Interval: +/- 0.485
    !!         T-Statistic:  5.158
    !!         P-Value:  0.000
    !! Coefficient 2:
    !!         Standard Error:  1.021
    !!         Confidence Interval: +/- 2.154
    !!         T-Statistic:  0.437
    !!         P-Value:  0.667
    !! Coefficient 3:
    !!         Standard Error:  1.204
    !!         Confidence Interval: +/- 2.540
    !!         T-Statistic: -0.102
    !!         P-Value:  0.920
    !! Coefficient 4:
    !!         Standard Error:  0.395
    !!         Confidence Interval: +/- 0.834
    !!         T-Statistic:  2.694
    !!         P-Value:  0.015
    !! Convergence Information:
    !!         Iteration Count: 9
    !!         Function Evaluation Count: 47
    !!         Converge on Gradient: F (0.215E-07)
    !!         Converge on Solution Change: T (0.466E-07)
    !!         Converge on Parameter Change: F (0.973E-01)
    !!
    !! Linear Model:
    !!         params(1):  1.18661
    !!         params(2):  0.44661
    !!         params(3): -0.12232
    !!         params(4):  1.06476
    !! Statistics:
    !! Coefficient 1:
    !!         Standard Error:  0.230
    !!         Confidence Interval: +/- 0.485
    !!         T-Statistic:  5.158
    !!         P-Value:  0.000
    !! Coefficient 2:
    !!         Standard Error:  1.021
    !!         Confidence Interval: +/- 2.154
    !!         T-Statistic:  0.437
    !!         P-Value:  0.667
    !! Coefficient 3:
    !!         Standard Error:  1.204
    !!         Confidence Interval: +/- 2.540
    !!         T-Statistic: -0.102
    !!         P-Value:  0.920
    !! Coefficient 4:
    !!         Standard Error:  0.395
    !!         Confidence Interval: +/- 0.834
    !!         T-Statistic:  2.694
    !!         P-Value:  0.015
    !! @endcode
    interface nonlinear_least_squares
        module procedure :: nonlinear_least_squares_1
    end interface

! ------------------------------------------------------------------------------
end module
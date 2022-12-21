! fstats.f90

module fstats
    use iso_fortran_env
    implicit none
    private
    public :: beta
    public :: regularized_beta
    public :: incomplete_beta
    public :: digamma
    public :: incomplete_gamma_upper
    public :: incomplete_gamma_lower
    public :: distribution
    public :: distribution_function
    public :: distribution_property
    public :: t_distribution
    public :: normal_distribution
    public :: f_distribution


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

! ------------------------------------------------------------------------------
end module
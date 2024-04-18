module fstats_distributions
    use iso_fortran_env
    use ieee_arithmetic
    use fstats_special_functions
    use fstats_helper_routines
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

    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

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
        procedure, public :: area => dist_area
            !! Computes the area under the PDF curve up to the value specified.
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

! ------------------------------------------------------------------------------
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
! ------------------------------------------------------------------------------
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

! ------------------------------------------------------------------------------
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

! ------------------------------------------------------------------------------
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

! ------------------------------------------------------------------------------
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

contains
! ------------------------------------------------------------------------------
pure elemental function dist_area(this, x) result(rst)
    !! Computes the area under the PDF curve up to the value of X specified.
    class(distribution), intent(in) :: this
        !! The distribution object.
    real(real64), intent(in) :: x
        !! The upper parameter limit.
    real(real64) :: rst
        !! The requested area.

    ! Local Variables
    integer(int32), parameter :: maxiter = 100
    real(real64), parameter :: tol = 1.0d-6
    integer(int32) :: i
    real(real64) :: f, df, h, twoh, dy

    ! Process
    !
    ! We use a simplified Newton's method to solve for the independent variable
    ! of the CDF function
    h = 1.0d-6
    twoh = 2.0d0 * h
    rst = 0.5d0 ! just an initial guess
    do i = 1, maxiter
        ! Compute the CDF and its derivative at y
        f = this%cdf(rst) - x
        df = (this%cdf(rst + h) - this%cdf(rst - h)) / twoh
        dy = f / df
        rst = rst - dy
        if (abs(dy) < tol) exit
    end do
end function

! ******************************************************************************
! STUDENT'S T-DISTRIBUTION
! ------------------------------------------------------------------------------
! REF: https://en.wikipedia.org/wiki/Student%27s_t-distribution
pure elemental function td_pdf(this, x) result(rst)
    !! Computes the probability density function.
    !!
    !! The PDF for Student's T-Distribution is given as 
    !! $$ f(t) = \frac{ \Gamma \left( \frac{\nu + 1}{2} \right) }
    !! { \sqrt{\nu \pi} \Gamma \left( \frac{\nu}{2} \right) } 
    !! \left( 1 + \frac{t^2}{\nu} \right)^{-(\nu + 1) / 2} $$.
    class(t_distribution), intent(in) :: this
        !! The t_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    ! Process
    rst = gamma((this%dof + 1.0d0) / 2.0d0) / &
        (sqrt(this%dof * pi) * gamma(this%dof / 2.0d0)) *&
        (1.0d0 + x**2 / this%dof)**(-0.5d0 * (1.0d0 + this%dof))
end function

! ------------------------------------------------------------------------------
pure elemental function td_cdf(this, x) result(rst)
    !! Computes the cumulative distribution function.
    !!
    !! The CDF for Student's T-Distribution is given as
    !! $$ F(t) = \int_{-\infty}^{t} f(u) \,du = 1 - \frac{1}{2} I_{x(t)}
    !! \left( \frac{\nu}{2}, \frac{1}{2} \right) $$
    !! where $$ x(t) = \frac{\nu}{\nu + t^2} $$.
    class(t_distribution), intent(in) :: this
        !! The t_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    ! Process
    real(real64) :: t
    t = this%dof / (this%dof + x**2)
    rst = 1.0d0 - 0.5d0 * regularized_beta(0.5d0 * this%dof, 0.5d0, t)
    if (x < 0) rst = 1.0d0 - rst
end function

! ------------------------------------------------------------------------------
pure function td_mean(this) result(rst)
    !! Computes the mean of the distribution.
    class(t_distribution), intent(in) :: this
        !! The t_distribution object.
    real(real64) :: rst
        !! The mean.

    ! Process
    if (this%dof < 1.0d0) then
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    else
        rst = 0.0d0
    end if
end function

! ------------------------------------------------------------------------------
pure function td_median(this) result(rst)
    !! Computes the median of the distribution.
    class(t_distribution), intent(in) :: this
        !! The t_distribution object.
    real(real64) :: rst

    ! Process
    rst = 0.0d0
end function

! ------------------------------------------------------------------------------
pure function td_mode(this) result(rst)
    !! Computes the mode of the distribution.
    class(t_distribution), intent(in) :: this
        !! The t_distribution object.
    real(real64) :: rst
        !! The mode.

    ! Process
    rst = 0.0d0
end function

! ------------------------------------------------------------------------------
pure function td_variance(this) result(rst)
    !! Computes the variance of the distribution.
    class(t_distribution), intent(in) :: this
        !! The t_distribution object.
    real(real64) :: rst
        !! The variance.

    ! Process
    if (this%dof <= 1.0d0) then
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    else if (this%dof > 1.0d0 .and. this%dof <= 2.0d0) then
        rst = ieee_value(rst, IEEE_POSITIVE_INF)
    else
        rst = this%dof / (this%dof - 2.0d0)
    end if
end function

! ******************************************************************************
! NORMAL DISTRIBUTION
! ------------------------------------------------------------------------------
pure elemental function nd_pdf(this, x) result(rst)
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

    rst = exp(-0.5d0 * ((x - this%mean_value) / this%standard_deviation)**2) / &
        (this%standard_deviation * sqrt(2.0d0 * pi))
end function

! ------------------------------------------------------------------------------
pure elemental function nd_cdf(this, x) result(rst)
    !! Computes the cumulative distribution function.
    !!
    !! The CDF for a normal distribution is given as 
    !! $$ F(x) = \frac{1}{2} \left( 1 + erf \left( \frac{x - \mu}
    !! {\sigma \sqrt{2}} \right) \right) $$.
    class(normal_distribution), intent(in) :: this
        !! The normal_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    rst = 0.5d0 * (1.0d0 + erf((x - this%mean_value) / &
        (this%standard_deviation * sqrt(2.0d0))))
end function

! ------------------------------------------------------------------------------
pure function nd_mean(this) result(rst)
    !! Computes the mean of the distribution.
    class(normal_distribution), intent(in) :: this
        !! The normal_distribution object.
    real(real64) :: rst
        !! The mean
    rst = this%mean_value
end function

! ------------------------------------------------------------------------------
pure function nd_median(this) result(rst)
    !! Computes the median of the distribution.
    class(normal_distribution), intent(in) :: this
        !! The normal_distribution object.
    real(real64) :: rst
        !! The median.
    rst = this%mean_value
end function

! ------------------------------------------------------------------------------
pure function nd_mode(this) result(rst)
    !! Computes the mode of the distribution.
    class(normal_distribution), intent(in) :: this
        !! The normal_distribution object.
    real(real64) :: rst
        !! The mode.
    rst = this%mean_value
end function

! ------------------------------------------------------------------------------
pure function nd_variance(this) result(rst)
    !! Computes the variance of the distribution.
    class(normal_distribution), intent(in) :: this
        !! The normal_distribution object.
    real(real64) :: rst
        !! The variance.
    rst = this%standard_deviation**2
end function

! ------------------------------------------------------------------------------
subroutine nd_standardize(this)
    !! Standardizes the normal distribution to a mean of 0 and a 
    !! standard deviation of 1.
    class(normal_distribution), intent(inout) :: this
        !! The normal_distribution object.
    this%mean_value = 0.0d0
    this%standard_deviation = 1.0d0
end subroutine

! ******************************************************************************
! F DISTRIBUTION
! ------------------------------------------------------------------------------
pure elemental function fd_pdf(this, x) result(rst)
    !! Computes the probability density function.
    !!
    !! The PDF for a F distribution is given as 
    !! $$ f(x) = 
    !! \sqrt{ \frac{ (d_1 x)^{d_1} d_{2}^{d_2} }{ (d_1 x + d_2)^{d_1 + d_2} } } 
    !! \frac{1}{x \beta \left( \frac{d_1}{2}, \frac{d_2}{2} \right) } $$.
    class(f_distribution), intent(in) :: this
        !! The f_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    ! Process
    real(real64) :: d1, d2
    d1 = this%d1
    d2 = this%d2
    rst = (1.0d0 / beta(0.5d0 * d1, 0.5d0 * d2)) * (d1 / d2)**(0.5d0 * d1) * &
        x**(0.5d0 * d1 - 1.0d0) * (1.0d0 + d1 * x/ d2)**(-0.5d0 * (d1 + d2))
end function

! ------------------------------------------------------------------------------
pure elemental function fd_cdf(this, x) result(rst)
    !! Computes the cumulative distribution function.
    !!
    !! The CDF for a F distribution is given as 
    !! $$ F(x) = I_{d_1 x/(d_1 x + d_2)} \left( \frac{d_1}{2}, 
    !! \frac{d_2}{2} \right) $$.
    class(f_distribution), intent(in) :: this
        !! The f_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    ! Process
    real(real64) :: d1, d2
    d1 = this%d1
    d2 = this%d2
    rst = regularized_beta(0.5d0 * d1, 0.5d0 * d2, d1 * x / (d1 * x + d2))
end function

! ------------------------------------------------------------------------------
pure function fd_mean(this) result(rst)
    !! Computes the mean of the distribution.
    class(f_distribution), intent(in) :: this
        !! The f_distribution object.
    real(real64) :: rst
        !! The mean.

    ! Process
    if (this%d2 > 2.0d0) then
        rst = this%d2 / (this%d2 - 2.0d0)
    else
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    end if
end function

! ------------------------------------------------------------------------------
pure function fd_median(this) result(rst)
    !! Computes the median of the distribution.
    class(f_distribution), intent(in) :: this
        !! The f_distribution object.
    real(real64) :: rst
        !! The median.
    rst = ieee_value(rst, IEEE_QUIET_NAN)
end function

! ------------------------------------------------------------------------------
pure function fd_mode(this) result(rst)
    !! Computes the mode of the distribution.
    class(f_distribution), intent(in) :: this
        !! The f_distribution object.
    real(real64) :: rst
        !! The mode.

    ! Process
    if (this%d1 > 2.0d0) then
        rst = ((this%d1 - 2.0d0) / this%d1) * (this%d2 / (this%d2 + 2.0d0))
    else
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    end if
end function

! ------------------------------------------------------------------------------
pure function fd_variance(this) result(rst)
    !! Computes the variance of the distribution.
    class(f_distribution), intent(in) :: this
        !! The f_distribution object.
    real(real64) :: rst
        !! The variance.

    ! Process
    real(real64) :: d1, d2
    d1 = this%d1
    d2 = this%d2
    if (d2 > 4.0d0) then
        rst = (2.0d0 * d2**2 * (d1 + d2 - 2.0d0)) / &
            (d1 * (d2 - 2.0d0)**2 * (d2 - 4.0d0))
    else
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    end if
end function

! ******************************************************************************
! CHI-SQUARED DISTRIBUTION
! ------------------------------------------------------------------------------
pure elemental function cs_pdf(this, x) result(rst)
    !! Computes the probability density function.
    !!
    !! The PDF for a Chi-squared distribution is given as 
    !! $$ f(x) = \frac{x^{k/2 - 1} \exp{-x / 2}} {2^{k / 2} 
    !! \Gamma \left( \frac{k}{2} \right)} $$.
    class(chi_squared_distribution), intent(in) :: this
        !! The chi_squared_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    ! Local Variables
    real(real64) :: arg

    ! Process
    arg = 0.5d0 * this%dof
    rst = 1.0d0 / (2.0d0**arg * gamma(arg)) * x**(arg - 1.0d0) * exp(-0.5d0 * x)
end function

! ------------------------------------------------------------------------------
pure elemental function cs_cdf(this, x) result(rst)
    !! Computes the cumulative distribution function.
    !!
    !! The CDF for a Chi-squared distribution is given as 
    !! $$ F(x) = \frac{ \gamma \left( \frac{k}{2}, \frac{x}{2} \right) }
    !! { \Gamma \left( \frac{k}{2} \right)} $$.
    class(chi_squared_distribution), intent(in) :: this
        !! The chi_squared_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    ! Local Variables
    real(real64) :: arg

    ! Process
    arg = 0.5d0 * this%dof
    rst = incomplete_gamma_lower(arg, 0.5d0 * x) / gamma(arg)
end function

! ------------------------------------------------------------------------------
pure function cs_mean(this) result(rst)
    !! Computes the mean of the distribution.
    class(chi_squared_distribution), intent(in) :: this
        !! The chi_squared_distribution object.
    real(real64) :: rst
        !! The mean.

    ! Process
    rst = real(this%dof, real64)
end function

! ------------------------------------------------------------------------------
pure function cs_median(this) result(rst)
    !! Computes the median of the distribution.
    class(chi_squared_distribution), intent(in) :: this
        !! The chi_squared_distribution object.
    real(real64) :: rst
        !! The median.

    ! Process
    rst = this%dof * (1.0d0 - 2.0d0 / (9.0d0 * this%dof))**3
end function

! ------------------------------------------------------------------------------
pure function cs_mode(this) result(rst)
    !! Computes the mode of the distribution.
    class(chi_squared_distribution), intent(in) :: this
        !! The chi_squared_distribution object.
    real(real64) :: rst
        !! The mode.

    ! Process
    rst = max(this%dof - 2.0d0, 0.0d0)
end function

! ------------------------------------------------------------------------------
pure function cs_variance(this) result(rst)
    !! Computes the variance of the distribution.
    class(chi_squared_distribution), intent(in) :: this
        !! The chi_squared_distribution object.
    real(real64) :: rst
        !! The variance.

    ! Process
    rst = 2.0d0 * this%dof
end function

! ******************************************************************************
! BINOMIAL DISTRIBUTION
! ------------------------------------------------------------------------------
pure elemental function bd_pdf(this, x) result(rst)
    !! Computes the probability mass function.
    !!
    !! The PMF for a binomial distribution is given as 
    !! $$ f(k,n,p) = \frac{n!}{k! \left( n - k! \right)} p^k 
    !! \left( 1 - p \right)^{n-k} $$.
    class(binomial_distribution), intent(in) :: this
        !! The binomial_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.  This parameter
        !! is the number k successes in the n independent trials.  As
        !! such, this parameter must exist on the set [0, n].
    real(real64) :: rst
        !! The value of the function.

    ! Local Variables
    real(real64) :: dn

    ! Process
    dn = real(this%n, real64)
    rst = (factorial(dn) / (factorial(x) * factorial(dn - x))) * (this%p**x) * (1.0d0 - this%p)**(dn - x)
end function

! ------------------------------------------------------------------------------
pure elemental function bd_cdf(this, x) result(rst)
    !! Computes the cumulative distribution funtion.
    !!
    !! The CDF for a binomial distribution is given as 
    !! $$ F(k,n,p) = I_{1-p} \left( n - k, 1 + k \right) $$, which is simply
    !! the regularized incomplete beta function.
    class(binomial_distribution), intent(in) :: this
        !! The binomial_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.  This parameter
        !! is the number k successes in the n independent trials.  As
        !! such, this parameter must exist on the set [0, n].
    real(real64) :: rst
        !! The value of the function.

    ! Local Variables
    real(real64) :: dn

    ! Process
    dn = real(this%n, real64)
    rst = regularized_beta(dn - x, x + 1.0d0, 1.0d0 - this%p)
end function

! ------------------------------------------------------------------------------
pure function bd_mean(this) result(rst)
    !! Computes the mean of the distribution.
    class(binomial_distribution), intent(in) :: this
        !! The binomial_distribution object.
    real(real64) :: rst
        !! The mean.

    rst = real(this%n * this%p, real64)
end function

! ------------------------------------------------------------------------------
pure function bd_median(this) result(rst)
    !! Computes the median of the distribution.
    class(binomial_distribution), intent(in) :: this
        !! The binomial_distribution object.
    real(real64) :: rst
        !! The median.

    rst = real(this%n * this%p, real64)
end function

! ------------------------------------------------------------------------------
pure function bd_mode(this) result(rst)
    !! Computes the mode of the distribution.
    class(binomial_distribution), intent(in) :: this
        !! The binomial_distribution object.
    real(real64) :: rst
        !! The mode.

    rst = (this%n + 1.0d0) * this%p
end function

! ------------------------------------------------------------------------------
pure function bd_variance(this) result(rst)
    !! Computes the variance of the distribution.
    class(binomial_distribution), intent(in) :: this
        !! The binomial_distribution object.
    real(real64) :: rst
        !! The variance.

    rst = this%n * this%p * (1.0d0 - this%p)
end function

! ------------------------------------------------------------------------------
end module
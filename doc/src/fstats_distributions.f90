module fstats_distributions
    use iso_fortran_env
    use ieee_arithmetic
    use fstats_special_functions
    use fstats_helper_routines
    use ferror
    use fstats_errors
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
    public :: multivariate_distribution
    public :: multivariate_distribution_function
    public :: multivariate_normal_distribution
    public :: log_normal_distribution

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
        procedure, public :: standardized_variable => dist_std_var
            !! Computes the standardized variable for the distribution.
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

! ------------------------------------------------------------------------------
    type, extends(distribution) :: log_normal_distribution
        !! Defines a normal distribution.
        real(real64) :: standard_deviation
            !! The standard deviation of the distribution.
        real(real64) :: mean_value
            !! The mean value of the distribution.
    contains
        procedure, public :: pdf => lnd_pdf
        procedure, public :: cdf => lnd_cdf
        procedure, public :: mean => lnd_mean
        procedure, public :: median => lnd_median
        procedure, public :: mode => lnd_mode
        procedure, public :: variance => lnd_variance
    end type

! ******************************************************************************
! MULTIVARIATE DISTRIBUTIONS
! ------------------------------------------------------------------------------
    type, abstract :: multivariate_distribution
        !! Defines a multivariate probability distribution.
    contains
        procedure(multivariate_distribution_function), deferred, pass :: pdf
            !! Computes the probability density function.
    end type

    interface
        pure function multivariate_distribution_function(this, x) result(rst)
            !! Defines an interface for a multivariate probability distribution
            !! function.
            use iso_fortran_env, only : real64
            import multivariate_distribution
            class(multivariate_distribution), intent(in) :: this
                !! The distribution object.
            real(real64), intent(in), dimension(:) :: x
                !! The values at which to evaluate the function.
            real(real64) :: rst
                !! The value of the function.
        end function
    end interface

! ------------------------------------------------------------------------------
    type, extends(multivariate_distribution) :: multivariate_normal_distribution
        !! Defines a multivariate normal (Gaussian) distribution.
        real(real64), private, allocatable, dimension(:) :: m_means
            !! An N-element array of mean values.
        real(real64), private, allocatable, dimension(:,:) :: m_cov
            !! The N-by-N covariance matrix.  This matrix must be 
            !! positive-definite.
        real(real64), private, allocatable, dimension(:,:) :: m_cholesky
            !! The N-by-N Cholesky factored form (lower) of the covariance
            !! matrix.
        real(real64), private, allocatable, dimension(:,:) :: m_covInv
            !! The N-by-N inverse of the covariance matrix.
        real(real64), private :: m_covDet
            !! The determinant of the covariance matrix.
    contains
        procedure, public :: initialize => mvnd_init
        procedure, public :: pdf => mvnd_pdf
        procedure, public :: get_means => mvnd_get_means
        procedure, public :: set_means => mvnd_update_mean
        procedure, public :: get_covariance => mvnd_get_covariance
        procedure, public :: get_cholesky_factored_matrix => mvnd_get_cholesky
    end type

contains
! ------------------------------------------------------------------------------
pure elemental function dist_std_var(this, x) result(rst)
    !! Computes the standardized variable for the distribution.
    class(distribution), intent(in) :: this
        !! The distribution object.
    real(real64), intent(in) :: x
        !! The value of interest.
    real(real64) :: rst
        !! The result.

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

! ******************************************************************************
! MULTIVARIATE NORMAL DISTRIBUTION
! ------------------------------------------------------------------------------
subroutine mvnd_init(this, mu, sigma, err)
    use linalg, only : cholesky_factor
    !! Initializes the multivariate normal distribution by defining the mean
    !! values and covariance matrix.
    class(multivariate_normal_distribution), intent(inout) :: this
        !! The multivariate_normal_distribution object.
    real(real64), intent(in), dimension(:) :: mu
        !! An N-element array containing the mean values.
    real(real64), intent(in), dimension(:,:) :: sigma
        !! The N-by-N covariance matrix.  The PDF exists only if this matrix
        !! is positive-definite; therefore, the positive-definite constraint 
        !! is checked within this routine and enforced.  An error is thrown if
        !! the supplied matrix is not positive-definite.
    class(errors), intent(inout), optional, target :: err
        !! The error handling object.

    ! Local Variables
    integer(int32) :: n, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(mu)

    ! Input Checking
    if (size(sigma, 1) /= n .or. size(sigma, 2) /= n) then
        call report_matrix_size_error(errmgr, "mvnd_init", "sigma", n, n, &
            size(sigma, 1), size(sigma, 2))
        return
    end if

    ! Store the matrices
    this%m_means = mu
    this%m_cov = sigma
    if (allocated(this%m_covInv)) then
        if (size(this%m_covInv, 1) /= n .or. size(this%m_covInv, 2) /= n) then
            deallocate(this%m_covInv)
            allocate(this%m_covInv(n, n), stat = flag)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_covInv(n, n), stat = flag)
        if (flag /= 0) go to 10
    end if
    if (allocated(this%m_cholesky)) then
        if (size(this%m_cholesky, 1) /= n .or. size(this%m_cholesky, 2) /= n) then
            deallocate(this%m_cholesky)
            allocate(this%m_cholesky(n, n), stat = flag)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_cholesky(n, n), stat = flag, source = sigma)
        if (flag /= 0) go to 10
    end if

    ! Compute the Cholesky factorization of the covariance matrix
    call cholesky_factor(this%m_cholesky, upper = .false., err = errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute the inverse and determinant
    call populate_identity(this%m_covInv)
    call cholesky_inverse(this%m_cholesky, this%m_covInv)
    this%m_covDet = cholesky_determinant(this%m_cholesky)

    ! End
    return

    ! Memory Error Handling
10  continue
    call report_memory_error(errmgr, "mvnd_init", flag)
    return
end subroutine

! ------------------------------------------------------------------------------
pure function mvnd_pdf(this, x) result(rst)
    !! Evaluates the PDF for the multivariate normal distribution.
    class(multivariate_normal_distribution), intent(in) :: this
        !! The multivariate_normal_distribution object.
    real(real64), intent(in), dimension(:) :: x
        !! The values at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    ! Local Variables
    integer(int32) :: n
    real(real64) :: arg
    real(real64), allocatable, dimension(:) :: delta, prod

    ! Process
    n = size(x)
    delta = x - this%m_means
    prod = matmul(this%m_covInv, delta) ! prod = inv(sigma) * (x - mu)
    arg = dot_product(delta, prod)      ! arg = (x - mu)**T * prod
    rst = exp(-0.5d0 * arg) / sqrt((2.0d0 * pi)**n * this%m_covDet)
end function

! ------------------------------------------------------------------------------
subroutine mvnd_update_mean(this, x, err)
    !! Updates the mean value array.
    class(multivariate_normal_distribution), intent(inout) :: this
        !! The multivariate_normal_distribution object.
    real(real64), intent(in), dimension(:) :: x
        !! The N-element array of new mean values.
    class(errors), intent(inout), optional, target :: err
        !! The error handling object.  This is referenced only in the event that
        !! the size of x is not compatible with the existing state.

    ! Local Variables
    integer(int32) :: n, nc, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)
    nc = size(this%m_means)

    ! Process
    if (.not.allocated(this%m_means)) then
        ! This is an initial set-up - just store the values and be done
        allocate(this%m_means(n), stat = flag, source = x)
        if (flag /= 0) then
            call report_memory_error(errmgr, "mvnd_update_mean", flag)
            return
        end if
        return
    end if

    ! Else, ensure the array is of the correct size before updating
    if (n /= nc) then
        call report_array_size_error(errmgr, "mvnd_update_mean", "x", nc, n)
        return
    end if
    this%m_means = x
end subroutine

! ------------------------------------------------------------------------------
pure function mvnd_get_means(this) result(rst)
    !! Gets the mean values of the distribution.
    class(multivariate_normal_distribution), intent(in) :: this
        !! The multivariate_normal_distribution object.
    real(real64), allocatable, dimension(:) :: rst
        !! The mean values.

    ! Process
    integer(int32) :: n
    if (allocated(this%m_means)) then
        n = size(this%m_means)
        allocate(rst(n), source = this%m_means)
    else
        allocate(rst(0))
    end if
end function

! ------------------------------------------------------------------------------
pure function mvnd_get_covariance(this) result(rst)
    !! Gets the covariance matrix of the distribution.
    class(multivariate_normal_distribution), intent(in) :: this
        !! The multivariate_normal_distribution object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The covariance matrix.

    ! Process
    integer(int32) :: n
    if (allocated(this%m_cov)) then
        n = size(this%m_cov, 1)
        allocate(rst(n, n), source = this%m_cov)
    else
        allocate(rst(0, 0))
    end if
end function

! ------------------------------------------------------------------------------
pure function mvnd_get_cholesky(this) result(rst)
    !! Gets the lower triangular form of the Cholesky factorization of the
    !! covariance matrix of the distribution.
    class(multivariate_normal_distribution), intent(in) :: this
        !! The multivariate_normal_distribution object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The Cholesky factored matrix.

    ! Process
    integer(int32) :: n
    if (allocated(this%m_cholesky)) then
        n = size(this%m_cholesky, 1)
        allocate(rst(n, n), source = this%m_cholesky)
    else
        allocate(rst(0, 0))
    end if
end function

! ******************************************************************************
! LOG NORMAL DISTRIBUTION
! ------------------------------------------------------------------------------
pure elemental function lnd_pdf(this, x) result(rst)
    !! Computes the probability density function.
    !!
    !! The PDF for a log-normal distribution is given as
    !! $$ f(x) = \frac{1}{x \sigma \sqrt{2 \pi}} \exp{\left(- \frac{\left( 
    !! \ln{x} - \mu \right)^2}{2 \sigma^2} \right)} $$
    class(log_normal_distribution), intent(in) :: this
        !! The log_normal_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    rst = exp(-(log(x) - this%mean_value)**2 / &
        (2.0d0 * this%standard_deviation**2)) / &
        (x * this%standard_deviation * sqrt(2.0d0 * pi))
end function

! ------------------------------------------------------------------------------
pure elemental function lnd_cdf(this, x) result(rst)
    !! Computes the cumulative distribution function.
    !!
    !! The CDF for a log-normal distribution is given as
    !! $$ F(x) = \frac{1}{2} \left(1 + erf\left( \frac{\ln{x} - \mu}
    !! {\sigma \sqrt{2}} \right) \right) $$
    class(log_normal_distribution), intent(in) :: this
        !! The log_normal_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    rst = 0.5d0 * (1.0d0 + erf((log(x) - this%standard_deviation) / &
        (this%standard_deviation * sqrt(2.0d0))))
end function

! ------------------------------------------------------------------------------
pure function lnd_mean(this) result(rst)
    !! Computes the mean of the distribution
    class(log_normal_distribution), intent(in) :: this
        !! The log_normal distribution object.
    real(real64) :: rst
        !! The mean
    rst = exp(this%mean_value + 0.5d0 * this%standard_deviation**2)
end function

! ------------------------------------------------------------------------------
pure function lnd_median(this) result(rst)
    !! Computes the median of the distribution
    class(log_normal_distribution), intent(in) :: this
        !! The log_normal distribution object.
    real(real64) :: rst
        !! The median
    rst = exp(this%mean_value)
end function

! ------------------------------------------------------------------------------
pure function lnd_mode(this) result(rst)
    !! Computes the mode of the distribution
    class(log_normal_distribution), intent(in) :: this
        !! The log_normal distribution object.
    real(real64) :: rst
        !! The mode
    rst = exp(this%mean_value - this%standard_deviation**2)
end function

! ------------------------------------------------------------------------------
pure function lnd_variance(this) result(rst)
    !! Computes the variance of the distribution
    class(log_normal_distribution), intent(in) :: this
        !! The log_normal distribution object.
    real(real64) :: rst
        !! The variance
    rst = (exp(this%standard_deviation**2) - 1.0d0) * &
        exp(2.0d0 * this%mean_value + this%standard_deviation**2)
end function

! ******************************************************************************
! SUPPORTING ROUTINES
! ------------------------------------------------------------------------------
subroutine cholesky_inverse(x, u)
    use linalg, only : solve_triangular_system
    !! Computes the inverse of a Cholesky-factored matrix.
    real(real64), intent(in), dimension(:,:) :: x
        !! The lower-triangular Cholesky factored matrix.
    real(real64), intent(inout), dimension(:,:) :: u
        !! On input, an N-by-N identity matrix.  On output, the N-by-N inverted
        !! matrix.

    ! To compute the inverse of a Cholesky factored matrix (L) consider the 
    ! following:
    !
    ! A = L * L**T
    !
    ! (L * L**T) * inv(A) = I, where I is an identity matrix
    !
    ! First, solve L * U = I, for the N-by-N matrix U
    !
    ! And then solve L' * inv(A) = U for inv(A)

    ! Solve L * U = I for U
    call solve_triangular_system(.true., .false., .false., .true., 1.0d0, x, u)

    ! Solve L**T * inv(A) = U for inv(A)
    call solve_triangular_system(.true., .false., .true., .true., 1.0d0, x, u)
end subroutine

! ------------------------------------------------------------------------------
pure function cholesky_determinant(x) result(rst)
    !! Computes the determinant of a Cholesky factored (lower) matrix.
    real(real64), intent(in), dimension(:,:) :: x
        !! The lower-triangular Cholesky-factored matrix.
    real(real64) :: rst
        !! The determinant.

    ! Local Variables
    integer(int32) :: i, ep, n
    real(real64) :: temp

    ! Initialization
    n = size(x, 1)
    rst = 0.0d0

    ! Compute the product of the squares of the diagonal
    temp = 1.0d0
    ep = 0
    do i = 1, n
        temp = (x(i,i))**2 * temp
        if (temp == 0.0d0) then
            rst = 0.0d0
            return
        end if

        do while (abs(temp) < 1.0d0)
            temp = 1.0d1 * temp
            ep = ep - 1
        end do

        do while (abs(temp) > 1.0d1)
            temp = 1.0d-1 * temp
            ep = ep + 1
        end do
    end do
    rst = temp * (1.0d1)**ep
end function

! ------------------------------------------------------------------------------
subroutine populate_identity(x)
    !! Populates the supplied matrix as an identity matrix.
    real(real64), intent(inout), dimension(:,:) :: x

    ! Local Variables
    integer(int32) :: i, m, n, mn

    ! Process
    m = size(x, 1)
    n = size(x, 2)
    mn = min(m, n)
    x = 0.0d0
    do i = 1, mn
        x(i,i) = 1.0d0
    end do
end subroutine

! ------------------------------------------------------------------------------
end module
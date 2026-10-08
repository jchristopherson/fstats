module fstats_chi_squared_distribution
    use iso_fortran_env
    use ieee_arithmetic
    use fstats_distributions, only : distribution, probability_log
    use fstats_special_functions
    implicit none
    private
    public :: chi_squared_distribution

    type, extends(distribution) :: chi_squared_distribution
        !! Defines a Chi-squared distribution.
        integer(int32) :: dof
            !! The number of degrees of freedom.
    contains
        procedure, public :: log_pdf => cs_log_pdf
        procedure, public :: survival => cs_survival
        procedure, public :: log_cdf => cs_log_cdf
        procedure, public :: log_survival => cs_log_survival
        procedure, public :: pdf => cs_pdf
        procedure, public :: cdf => cs_cdf
        procedure, public :: mean => cs_mean
        procedure, public :: median => cs_median
        procedure, public :: mode => cs_mode
        procedure, public :: variance => cs_variance
        procedure, public :: defined_range => cs_range
        procedure, public :: recenter => cs_recenter
    end type

contains

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
    rst = exp(this%log_pdf(x))
end function

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
    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (this%dof <= 0 .or. ieee_is_nan(x)) return
    if (x <= 0.0d0) then
        rst = 0.0d0
    else
        rst = regularized_gamma_lower(arg, 0.5d0 * x)
    end if
end function

pure function cs_mean(this) result(rst)
    !! Computes the mean of the distribution.
    class(chi_squared_distribution), intent(in) :: this
        !! The chi_squared_distribution object.
    real(real64) :: rst
        !! The mean.

    ! Process
    rst = real(this%dof, real64)
end function

pure function cs_median(this) result(rst)
    !! Computes the median of the distribution.
    class(chi_squared_distribution), intent(in) :: this
        !! The chi_squared_distribution object.
    real(real64) :: rst
        !! The median.

    ! Process
    rst = this%dof * (1.0d0 - 2.0d0 / (9.0d0 * this%dof))**3
end function

pure function cs_mode(this) result(rst)
    !! Computes the mode of the distribution.
    class(chi_squared_distribution), intent(in) :: this
        !! The chi_squared_distribution object.
    real(real64) :: rst
        !! The mode.

    ! Process
    rst = max(this%dof - 2.0d0, 0.0d0)
end function

pure function cs_variance(this) result(rst)
    !! Computes the variance of the distribution.
    class(chi_squared_distribution), intent(in) :: this
        !! The chi_squared_distribution object.
    real(real64) :: rst
        !! The variance.

    ! Process
    rst = 2.0d0 * this%dof
end function

pure function cs_range(this) result(rst)
    !! Gets the defined range for the distribution.
    class(chi_squared_distribution), intent(in) :: this
        !! The chi_squared_distribution object.
    real(real64), dimension(2) :: rst
        !! The defined range of the probability distributions [0, infinity).  As
        !! using a value of infinity may cause issue, this routine returns
        !! huge(0.0d0) instead.

    rst = [0.0d0, huge(0.0d0)]
end function

subroutine cs_recenter(this, x)
    !! Recenters the distribution about the supplied value.
    class(chi_squared_distribution), intent(inout) :: this
        !! The chi_squared_distribution object.
    real(real64), intent(in) :: x
        !! The value about which to recenter.

    ! Process
    this%dof = floor(x)
end subroutine

pure elemental function cs_log_pdf(this, x) result(rst)
    class(chi_squared_distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point; discrete masses require integer-valued x.
    real(real64) :: rst
        !! Log density; outside support gives -infinity, invalid inputs give NaN.
        !! An integrable density singularity at a support endpoint gives +infinity.
    real(real64) :: first

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (ieee_is_nan(x)) return
    if (this%dof <= 0) return
    rst = ieee_value(0.0d0, ieee_negative_inf)
    if (x < 0.0d0 .or. .not.ieee_is_finite(x)) return
    first = 0.5d0 * real(this%dof, real64)
    if (x == 0.0d0) then
        if (first < 1.0d0) rst = ieee_value(0.0d0, ieee_positive_inf)
        if (first == 1.0d0) rst = -log(2.0d0)
        return
    end if
    rst = (first - 1.0d0) * log(x) - 0.5d0 * x - first * log(2.0d0) - log_gamma(first)
end function

pure elemental function cs_survival(this, x) result(rst)
    class(chi_squared_distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point; discrete tails include masses strictly above x.
    real(real64) :: rst
        !! Probability on [0, 1]; invalid parameters or NaN x give NaN.

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (ieee_is_nan(x)) return
    if (this%dof <= 0) return
    rst = 1.0d0
    if (x <= 0.0d0) return
    rst = regularized_gamma_upper(0.5d0 * real(this%dof, real64), 0.5d0 * x)
end function

pure elemental function cs_log_cdf(this, x) result(rst)
    class(chi_squared_distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point.
    real(real64) :: rst
        !! Log probability; zero probability gives -infinity, invalid inputs NaN.

    rst = probability_log(this%cdf(x))
    if (ieee_is_nan(rst)) return
    if (x > 0.0d0) rst = log_regularized_gamma_lower(0.5d0 * real(this%dof, real64), 0.5d0 * x)
end function

pure elemental function cs_log_survival(this, x) result(rst)
    class(chi_squared_distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point.
    real(real64) :: rst
        !! Log probability; zero probability gives -infinity, invalid inputs NaN.

    rst = probability_log(this%survival(x))
    if (ieee_is_nan(rst)) return
    if (x > 0.0d0) rst = log_regularized_gamma_upper(0.5d0 * real(this%dof, real64), 0.5d0 * x)
end function

end module fstats_chi_squared_distribution

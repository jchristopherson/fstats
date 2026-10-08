module fstats_normal_distribution
    use iso_fortran_env
    use ieee_arithmetic
    use fstats_distributions, only : distribution, probability_log, pi => distribution_pi
    use fstats_special_functions
    implicit none
    private
    public :: normal_distribution

    type, extends(distribution) :: normal_distribution
        !! Defines a normal distribution.
        real(real64) :: standard_deviation
            !! The standard deviation of the distribution.
        real(real64) :: mean_value
            !! The mean value of the distribution.
    contains
        procedure, public :: pdf => nd_pdf
        procedure, public :: cdf => nd_cdf
        procedure, public :: log_pdf => nd_log_pdf
        procedure, public :: survival => nd_survival
        procedure, public :: log_cdf => nd_log_cdf
        procedure, public :: log_survival => nd_log_survival
        procedure, public :: mean => nd_mean
        procedure, public :: median => nd_median
        procedure, public :: mode => nd_mode
        procedure, public :: variance => nd_variance
        procedure, public :: standardize => nd_standardize
        procedure, public :: recenter => nd_recenter
    end type

contains

pure elemental function nd_log_pdf(this, x) result(rst)
    !! Computes a normal log density directly, without first evaluating its PDF.
    class(normal_distribution), intent(in) :: this
        !! Normal law with finite mean and finite, positive standard deviation.
    real(real64), intent(in) :: x
        !! Evaluation point; infinite x gives -infinity.
    real(real64) :: rst
        !! Log density; invalid parameters or NaN x give NaN.
    real(real64) :: standardized

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (ieee_is_nan(x)) return
    if (.not.ieee_is_finite(this%mean_value)) return
    if (.not.ieee_is_finite(this%standard_deviation)) return
    if (this%standard_deviation <= 0.0d0) return
    standardized = normal_standardized(x, this%mean_value, this%standard_deviation) / sqrt(2.0d0)
    if (abs(standardized) > sqrt(huge(1.0d0))) then
        rst = ieee_value(0.0d0, ieee_negative_inf)
    else
        rst = -standardized**2 - log(this%standard_deviation) - 0.5d0 * log(2.0d0 * pi)
    end if
end function

pure elemental function nd_survival(this, x) result(rst)
    !! Computes P(X > x) for a normal law using the complementary error function.
    class(normal_distribution), intent(in) :: this
        !! Normal law with finite mean and finite, positive standard deviation.
    real(real64), intent(in) :: x
        !! Evaluation point.
    real(real64) :: rst
        !! Upper-tail probability; invalid parameters or NaN x give NaN.

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (ieee_is_nan(x)) return
    if (.not.ieee_is_finite(this%mean_value)) return
    if (.not.ieee_is_finite(this%standard_deviation)) return
    if (this%standard_deviation <= 0.0d0) return
    rst = 0.5d0 * erfc(normal_standardized(x, this%mean_value, this%standard_deviation) / sqrt(2.0d0))
end function

pure elemental function nd_log_cdf(this, x) result(rst)
    !! Computes the normal log lower tail, including extreme negative tails.
    class(normal_distribution), intent(in) :: this
        !! Normal law with finite mean and finite, positive standard deviation.
    real(real64), intent(in) :: x
        !! Evaluation point.
    real(real64) :: rst
        !! Log CDF; invalid parameters or NaN x give NaN.

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (ieee_is_nan(x)) return
    if (.not.ieee_is_finite(this%mean_value)) return
    if (.not.ieee_is_finite(this%standard_deviation)) return
    if (this%standard_deviation <= 0.0d0) return
    rst = normal_log_tail(-normal_standardized(x, this%mean_value, this%standard_deviation))
end function

pure elemental function nd_log_survival(this, x) result(rst)
    !! Computes the normal log upper tail, including extreme positive tails.
    class(normal_distribution), intent(in) :: this
        !! Normal law with finite mean and finite, positive standard deviation.
    real(real64), intent(in) :: x
        !! Evaluation point.
    real(real64) :: rst
        !! Log survival probability; invalid parameters or NaN x give NaN.

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (ieee_is_nan(x)) return
    if (.not.ieee_is_finite(this%mean_value)) return
    if (.not.ieee_is_finite(this%standard_deviation)) return
    if (this%standard_deviation <= 0.0d0) return
    rst = normal_log_tail(normal_standardized(x, this%mean_value, this%standard_deviation))
end function

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

    rst = exp(this%log_pdf(x))
end function

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

    rst = exp(this%log_cdf(x))
end function

pure function nd_mean(this) result(rst)
    !! Computes the mean of the distribution.
    class(normal_distribution), intent(in) :: this
        !! The normal_distribution object.
    real(real64) :: rst
        !! The mean
    rst = this%mean_value
end function

pure function nd_median(this) result(rst)
    !! Computes the median of the distribution.
    class(normal_distribution), intent(in) :: this
        !! The normal_distribution object.
    real(real64) :: rst
        !! The median.
    rst = this%mean_value
end function

pure function nd_mode(this) result(rst)
    !! Computes the mode of the distribution.
    class(normal_distribution), intent(in) :: this
        !! The normal_distribution object.
    real(real64) :: rst
        !! The mode.
    rst = this%mean_value
end function

pure function nd_variance(this) result(rst)
    !! Computes the variance of the distribution.
    class(normal_distribution), intent(in) :: this
        !! The normal_distribution object.
    real(real64) :: rst
        !! The variance.
    rst = this%standard_deviation**2
end function

subroutine nd_standardize(this)
    !! Standardizes the normal distribution to a mean of 0 and a 
    !! standard deviation of 1.
    class(normal_distribution), intent(inout) :: this
        !! The normal_distribution object.
    this%mean_value = 0.0d0
    this%standard_deviation = 1.0d0
end subroutine

subroutine nd_recenter(this, x)
    !! Recenters the distribution about the supplied value.
    class(normal_distribution), intent(inout) :: this
        !! The normal_distribution object.
    real(real64), intent(in) :: x
        !! The value about which to recenter.

    ! Process
    this%mean_value = x
end subroutine

pure elemental function normal_standardized(x, mu, sigma) result(rst)
    !! Standardizes a normal evaluation point without an overflowing subtraction.
    real(real64), intent(in) :: x
        !! Evaluation point; NaN inputs propagate.
    real(real64), intent(in) :: mu
        !! Finite normal mean.
    real(real64), intent(in) :: sigma
        !! Finite, strictly positive normal standard deviation.
    real(real64) :: rst
        !! Standardized point (x - mu) / sigma.

    if ((x >= 0.0d0 .and. mu >= 0.0d0) .or. (x <= 0.0d0 .and. mu <= 0.0d0)) then
        rst = (x - mu) / sigma
    else
        rst = x / sigma - mu / sigma
    end if
end function

pure elemental function normal_log_tail(z) result(rst)
    !! Computes the standard-normal log upper tail.
    !! Uses ERFC centrally and an asymptotic Mills-ratio expansion for z >= 26,
    !! retaining finite log probabilities after the ordinary tail underflows.
    real(real64), intent(in) :: z
        !! Standardized evaluation point.
    real(real64) :: rst
        !! Log upper-tail probability; +infinity gives -infinity, NaN propagates.
    real(real64) :: term, total, previous, inverse_square
    integer(int32) :: iteration

    if (ieee_is_nan(z)) then
        rst = z
    else if (z < 0.0d0) then
        rst = stable_log1p(-0.5d0 * erfc(-z / sqrt(2.0d0)))
    else if (z < 26.0d0) then
        rst = probability_log(0.5d0 * erfc(z / sqrt(2.0d0)))
    else if (.not.ieee_is_finite(z) .or. z / sqrt(2.0d0) > sqrt(huge(1.0d0))) then
        rst = ieee_value(0.0d0, ieee_negative_inf)
    else
        inverse_square = (1.0d0 / z)**2
        term = 1.0d0
        total = term
        do iteration = 1, 100
            previous = term
            term = -term * real(2 * iteration - 1, real64) * inverse_square
            if (abs(term) >= abs(previous)) exit
            total = total + term
            if (abs(term) <= epsilon(total) * abs(total)) exit
        end do
        rst = -(z / sqrt(2.0d0))**2 - log(z) - 0.5d0 * log(2.0d0 * pi) + log(total)
    end if
end function

end module fstats_normal_distribution

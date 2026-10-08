module fstats_binomial_distribution
    use iso_fortran_env
    use ieee_arithmetic
    use fstats_distributions, only : distribution, probability_log
    use fstats_special_functions
    implicit none
    private
    public :: binomial_distribution

    type, extends(distribution) :: binomial_distribution
        !! Defines a binomial distribution.  The binomial distribution describes
        !! the probability p of getting k successes in n independent trials.
        integer(int32) :: n
            !! The number of independent trials.
        real(real64) :: p
            !! The success probability for each trial.  This parameter must
            !! exist on the set [0, 1].
    contains
        procedure, public :: log_pdf => bd_log_pdf
        procedure, public :: survival => bd_survival
        procedure, public :: log_cdf => bd_log_cdf
        procedure, public :: log_survival => bd_log_survival
        procedure, public :: pdf => bd_pdf
        procedure, public :: cdf => bd_cdf
        procedure, public :: mean => bd_mean
        procedure, public :: median => bd_median
        procedure, public :: mode => bd_mode
        procedure, public :: variance => bd_variance
        procedure, public :: defined_range => bd_range
        procedure, public :: recenter => bd_recenter
    end type

contains

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
    rst = exp(this%log_pdf(x))
end function

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
    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (this%n < 0 .or. .not.ieee_is_finite(this%p) .or. ieee_is_nan(x)) return
    if (this%p < 0.0d0 .or. this%p > 1.0d0) return
    if (x < 0.0d0) then
        rst = 0.0d0
    else if (x >= dn) then
        rst = 1.0d0
    else
        rst = regularized_beta(dn - aint(x), aint(x) + 1.0d0, 1.0d0 - this%p)
    end if
end function

pure function bd_mean(this) result(rst)
    !! Computes the mean of the distribution.
    class(binomial_distribution), intent(in) :: this
        !! The binomial_distribution object.
    real(real64) :: rst
        !! The mean.

    rst = real(this%n * this%p, real64)
end function

pure function bd_median(this) result(rst)
    !! Computes the median of the distribution.
    class(binomial_distribution), intent(in) :: this
        !! The binomial_distribution object.
    real(real64) :: rst
        !! The median.

    rst = real(this%n * this%p, real64)
end function

pure function bd_mode(this) result(rst)
    !! Computes the mode of the distribution.
    class(binomial_distribution), intent(in) :: this
        !! The binomial_distribution object.
    real(real64) :: rst
        !! The mode.

    rst = (this%n + 1.0d0) * this%p
end function

pure function bd_variance(this) result(rst)
    !! Computes the variance of the distribution.
    class(binomial_distribution), intent(in) :: this
        !! The binomial_distribution object.
    real(real64) :: rst
        !! The variance.

    rst = this%n * this%p * (1.0d0 - this%p)
end function

pure function bd_range(this) result(rst)
    !! Gets the defined range for the distribution.
    class(binomial_distribution), intent(in) :: this
        !! The binomial_distribution object.
    real(real64), dimension(2) :: rst
        !! The defined range of the probability distributions [0, infinity).  As
        !! using a value of infinity may cause issue, this routine returns
        !! huge(0.0d0) instead.

    rst = [0.0d0, huge(0.0d0)]
end function

subroutine bd_recenter(this, x)
    !! Recenters the distribution about the supplied value.
    class(binomial_distribution), intent(inout) :: this
        !! The binomial_distribution object.
    real(real64), intent(in) :: x
        !! The value about which to recenter.

    ! Process
    this%p = x / this%n
end subroutine

pure elemental function bd_log_pdf(this, x) result(rst)
    class(binomial_distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point; discrete masses require integer-valued x.
    real(real64) :: rst
        !! Log density; outside support gives -infinity, invalid inputs give NaN.
        !! An integrable density singularity at a support endpoint gives +infinity.
    real(real64) :: count

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (ieee_is_nan(x)) return
    if (this%n < 0 .or. .not.ieee_is_finite(this%p)) return
    if (this%p < 0.0d0 .or. this%p > 1.0d0) return
    rst = ieee_value(0.0d0, ieee_negative_inf)
    count = real(this%n, real64)
    if (x < 0.0d0 .or. x > count) return
    if (x /= aint(x)) return
    if (this%p == 0.0d0 .or. this%p == 1.0d0) then
        if (x == count * this%p) rst = 0.0d0
        return
    end if
    rst = log_gamma(count + 1.0d0) - log_gamma(x + 1.0d0) - log_gamma(count - x + 1.0d0) &
        + x * log(this%p) + (count - x) * stable_log1p(-this%p)
end function

pure elemental function bd_survival(this, x) result(rst)
    class(binomial_distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point; discrete tails include masses strictly above x.
    real(real64) :: rst
        !! Probability on [0, 1]; invalid parameters or NaN x give NaN.
    real(real64) :: count

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (ieee_is_nan(x)) return
    if (this%n < 0 .or. .not.ieee_is_finite(this%p)) return
    if (this%p < 0.0d0 .or. this%p > 1.0d0) return
    count = real(this%n, real64)
    if (x < 0.0d0) then
        rst = 1.0d0
    else if (x >= count) then
        rst = 0.0d0
    else
        rst = regularized_beta(aint(x) + 1.0d0, count - aint(x), this%p)
    end if
end function

pure elemental function bd_log_cdf(this, x) result(rst)
    class(binomial_distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point.
    real(real64) :: rst
        !! Log probability; zero probability gives -infinity, invalid inputs NaN.

    rst = probability_log(this%cdf(x))
    if (ieee_is_nan(rst)) return
    if (x >= 0.0d0 .and. x < real(this%n, real64)) &
        rst = log_regularized_beta(real(this%n, real64) - aint(x), aint(x) + 1.0d0, 1.0d0 - this%p)
end function

pure elemental function bd_log_survival(this, x) result(rst)
    class(binomial_distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point.
    real(real64) :: rst
        !! Log probability; zero probability gives -infinity, invalid inputs NaN.

    rst = probability_log(this%survival(x))
    if (ieee_is_nan(rst)) return
    if (x >= 0.0d0 .and. x < real(this%n, real64)) &
        rst = log_regularized_beta(aint(x) + 1.0d0, real(this%n, real64) - aint(x), this%p)
end function

end module fstats_binomial_distribution

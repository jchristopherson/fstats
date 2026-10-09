! SPDX-FileCopyrightText: 2022-2026 Jason Christopherson
! SPDX-License-Identifier: MIT
!
! MIT License
!
! Copyright (c) 2022-2026 Jason Christopherson
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module fstats_log_normal_distribution
    use iso_fortran_env
    use ieee_arithmetic
    use fstats_distributions, only : distribution, probability_log
    use fstats_special_functions
    use fstats_normal_distribution, only : normal_distribution
    implicit none
    private
    public :: log_normal_distribution

    type, extends(distribution) :: log_normal_distribution
        !! Defines a normal distribution.
        real(real64) :: standard_deviation
            !! The standard deviation of the distribution.
        real(real64) :: mean_value
            !! The mean value of the distribution.
    contains
        procedure, public :: pdf => lnd_pdf
        procedure, public :: cdf => lnd_cdf
        procedure, public :: log_pdf => lnd_log_pdf
        procedure, public :: survival => lnd_survival
        procedure, public :: log_cdf => lnd_log_cdf
        procedure, public :: log_survival => lnd_log_survival
        procedure, public :: mean => lnd_mean
        procedure, public :: median => lnd_median
        procedure, public :: mode => lnd_mode
        procedure, public :: variance => lnd_variance
        procedure, public :: defined_range => lnd_range
        procedure, public :: recenter => lnd_recenter
    end type

contains

pure elemental function lnd_log_pdf(this, x) result(rst)
    !! Computes a log-normal log density from the underlying normal log density.
    class(log_normal_distribution), intent(in) :: this
        !! Law with finite log-mean and finite, positive log-standard deviation.
    real(real64), intent(in) :: x
        !! Evaluation point in observation space, not log space.
    real(real64) :: rst
        !! Log density; x <= 0 gives -infinity, invalid parameters or NaN x NaN.
    type(normal_distribution) :: normal

    normal%mean_value = this%mean_value
    normal%standard_deviation = this%standard_deviation
    rst = normal%log_pdf(0.0d0)
    if (ieee_is_nan(rst) .or. ieee_is_nan(x)) then
        rst = ieee_value(0.0d0, ieee_quiet_nan)
        return
    end if
    if (x <= 0.0d0 .or. .not.ieee_is_finite(x)) then
        rst = ieee_value(0.0d0, ieee_negative_inf)
    else
        rst = normal%log_pdf(log(x)) - log(x)
    end if
end function

pure elemental function lnd_survival(this, x) result(rst)
    !! Computes the log-normal upper-tail probability P(X > x).
    class(log_normal_distribution), intent(in) :: this
        !! Law with finite log-mean and finite, positive log-standard deviation.
    real(real64), intent(in) :: x
        !! Evaluation point in observation space.
    real(real64) :: rst
        !! Survival probability; x <= 0 gives one, invalid inputs give NaN.

    rst = exp(this%log_survival(x))
end function

pure elemental function lnd_log_cdf(this, x) result(rst)
    !! Computes the log-normal log lower tail using the normal log-CDF.
    class(log_normal_distribution), intent(in) :: this
        !! Law with finite log-mean and finite, positive log-standard deviation.
    real(real64), intent(in) :: x
        !! Evaluation point in observation space.
    real(real64) :: rst
        !! Log CDF; x <= 0 gives -infinity, invalid inputs give NaN.
    type(normal_distribution) :: normal

    normal%mean_value = this%mean_value
    normal%standard_deviation = this%standard_deviation
    rst = normal%log_cdf(0.0d0)
    if (ieee_is_nan(rst) .or. ieee_is_nan(x)) then
        rst = ieee_value(0.0d0, ieee_quiet_nan)
    else if (x <= 0.0d0) then
        rst = ieee_value(0.0d0, ieee_negative_inf)
    else
        rst = normal%log_cdf(log(x))
    end if
end function

pure elemental function lnd_log_survival(this, x) result(rst)
    !! Computes the log-normal log upper tail using the normal log-survival.
    class(log_normal_distribution), intent(in) :: this
        !! Law with finite log-mean and finite, positive log-standard deviation.
    real(real64), intent(in) :: x
        !! Evaluation point in observation space.
    real(real64) :: rst
        !! Log survival probability; x <= 0 gives zero, invalid inputs give NaN.
    type(normal_distribution) :: normal

    normal%mean_value = this%mean_value
    normal%standard_deviation = this%standard_deviation
    rst = normal%log_survival(0.0d0)
    if (ieee_is_nan(rst) .or. ieee_is_nan(x)) then
        rst = ieee_value(0.0d0, ieee_quiet_nan)
    else if (x <= 0.0d0) then
        rst = 0.0d0
    else
        rst = normal%log_survival(log(x))
    end if
end function

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

    rst = exp(this%log_pdf(x))
end function

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

    rst = exp(this%log_cdf(x))
end function

pure function lnd_mean(this) result(rst)
    !! Computes the mean of the distribution
    class(log_normal_distribution), intent(in) :: this
        !! The log_normal distribution object.
    real(real64) :: rst
        !! The mean
    rst = exp(this%mean_value + 0.5d0 * this%standard_deviation**2)
end function

pure function lnd_median(this) result(rst)
    !! Computes the median of the distribution
    class(log_normal_distribution), intent(in) :: this
        !! The log_normal distribution object.
    real(real64) :: rst
        !! The median
    rst = exp(this%mean_value)
end function

pure function lnd_mode(this) result(rst)
    !! Computes the mode of the distribution
    class(log_normal_distribution), intent(in) :: this
        !! The log_normal distribution object.
    real(real64) :: rst
        !! The mode
    rst = exp(this%mean_value - this%standard_deviation**2)
end function

pure function lnd_variance(this) result(rst)
    !! Computes the variance of the distribution
    class(log_normal_distribution), intent(in) :: this
        !! The log_normal distribution object.
    real(real64) :: rst
        !! The variance
    rst = stable_expm1(this%standard_deviation**2) * &
        exp(2.0d0 * this%mean_value + this%standard_deviation**2)
end function

pure function lnd_range(this) result(rst)
    !! Gets the defined range for the distribution.
    class(log_normal_distribution), intent(in) :: this
        !! The log_normal_distribution object.
    real(real64), dimension(2) :: rst
        !! The defined range of the probability distributions [0, infinity).  As
        !! using a value of infinity may cause issue, this routine returns
        !! huge(0.0d0) instead.

    rst = [0.0d0, huge(0.0d0)]
end function

subroutine lnd_recenter(this, x)
    !! Recenters the distribution about the supplied value.
    class(log_normal_distribution), intent(inout) :: this
        !! The log_normal_distribution object.
    real(real64), intent(in) :: x
        !! The value about which to recenter.

    this%mean_value = x
end subroutine

end module fstats_log_normal_distribution

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

module fstats_f_distribution
    use iso_fortran_env
    use ieee_arithmetic
    use fstats_distributions, only : distribution, probability_log, softplus
    use fstats_special_functions
    implicit none
    private
    public :: f_distribution

    type, extends(distribution) :: f_distribution
        !! Defines an F-distribution.
        real(real64) :: d1
            !! The measure of degrees of freedom for the first data set.
        real(real64) :: d2
            !! The measure of degrees of freedom for the second data set.
    contains
        procedure, public :: log_pdf => fd_log_pdf
        procedure, public :: survival => fd_survival
        procedure, public :: log_cdf => fd_log_cdf
        procedure, public :: log_survival => fd_log_survival
        procedure, public :: pdf => fd_pdf
        procedure, public :: cdf => fd_cdf
        procedure, public :: mean => fd_mean
        procedure, public :: median => fd_median
        procedure, public :: mode => fd_mode
        procedure, public :: variance => fd_variance
        procedure, public :: defined_range => fd_range
        procedure, public :: recenter => fd_recenter
    end type

contains

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
    rst = exp(this%log_pdf(x))
end function

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
    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (.not.ieee_is_finite(d1) .or. .not.ieee_is_finite(d2)) return
    if (d1 <= 0.0d0 .or. d2 <= 0.0d0 .or. ieee_is_nan(x)) return
    if (x <= 0.0d0) then
        rst = 0.0d0
    else
        rst = regularized_beta(0.5d0 * d1, 0.5d0 * d2, &
            exp(-softplus(log(d2) - log(d1) - log(x))))
    end if
end function

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

pure function fd_median(this) result(rst)
    !! Computes the median of the distribution.
    class(f_distribution), intent(in) :: this
        !! The f_distribution object.
    real(real64) :: rst
        !! The median.
    rst = ieee_value(rst, IEEE_QUIET_NAN)
end function

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

pure function fd_range(this) result(rst)
    !! Gets the defined range for the distribution.
    class(f_distribution), intent(in) :: this
        !! The f_distribution object.
    real(real64), dimension(2) :: rst
        !! The defined range of the probability distributions [0, infinity).  As
        !! using a value of infinity may cause issue, this routine returns
        !! huge(0.0d0) instead.

    rst = [0.0d0, huge(0.0d0)]
end function

subroutine fd_recenter(this, x)
    !! Recenters the distribution about the supplied value.
    class(f_distribution), intent(inout) :: this
        !! The f_distribution object.
    real(real64), intent(in) :: x
        !! The value about which to recenter.

    ! Process
    this%d2 = 2.0d0 * x / (x - 1.0d0)
end subroutine

pure elemental function fd_log_pdf(this, x) result(rst)
    class(f_distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point; discrete masses require integer-valued x.
    real(real64) :: rst
        !! Log density; outside support gives -infinity, invalid inputs give NaN.
        !! An integrable density singularity at a support endpoint gives +infinity.
    real(real64) :: first, second, argument

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (ieee_is_nan(x)) return
    if (.not.ieee_is_finite(this%d1) .or. .not.ieee_is_finite(this%d2)) return
    if (this%d1 <= 0.0d0 .or. this%d2 <= 0.0d0) return
    rst = ieee_value(0.0d0, ieee_negative_inf)
    if (x < 0.0d0 .or. .not.ieee_is_finite(x)) return
    first = 0.5d0 * this%d1
    second = 0.5d0 * this%d2
    if (x == 0.0d0) then
        if (first < 1.0d0) rst = ieee_value(0.0d0, ieee_positive_inf)
        if (first == 1.0d0) rst = 0.0d0
        return
    end if
    argument = log(this%d1) - log(this%d2)
    rst = first * argument + (first - 1.0d0) * log(x) &
        - (first + second) * softplus(argument + log(x)) &
        - log_gamma(first) - log_gamma(second) + log_gamma(first + second)
end function

pure elemental function fd_survival(this, x) result(rst)
    class(f_distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point; discrete tails include masses strictly above x.
    real(real64) :: rst
        !! Probability on [0, 1]; invalid parameters or NaN x give NaN.
    real(real64) :: argument

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (ieee_is_nan(x)) return
    if (.not.ieee_is_finite(this%d1) .or. .not.ieee_is_finite(this%d2)) return
    if (this%d1 <= 0.0d0 .or. this%d2 <= 0.0d0) return
    rst = 1.0d0
    if (x <= 0.0d0) return
    argument = exp(-softplus(log(this%d1) + log(x) - log(this%d2)))
    rst = regularized_beta(0.5d0 * this%d2, 0.5d0 * this%d1, argument)
end function

pure elemental function fd_log_cdf(this, x) result(rst)
    class(f_distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point.
    real(real64) :: rst
        !! Log probability; zero probability gives -infinity, invalid inputs NaN.
    real(real64) :: argument

    rst = probability_log(this%cdf(x))
    if (ieee_is_nan(rst)) return
    if (x > 0.0d0) then
        argument = exp(-softplus(log(this%d2) - log(this%d1) - log(x)))
        rst = log_regularized_beta(0.5d0 * this%d1, 0.5d0 * this%d2, argument)
    end if
end function

pure elemental function fd_log_survival(this, x) result(rst)
    class(f_distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point.
    real(real64) :: rst
        !! Log probability; zero probability gives -infinity, invalid inputs NaN.
    real(real64) :: argument

    rst = probability_log(this%survival(x))
    if (ieee_is_nan(rst)) return
    if (x > 0.0d0) then
        argument = exp(-softplus(log(this%d1) + log(x) - log(this%d2)))
        rst = log_regularized_beta(0.5d0 * this%d2, 0.5d0 * this%d1, argument)
    end if
end function

end module fstats_f_distribution

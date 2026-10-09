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

module fstats_t_distribution
    use iso_fortran_env
    use ieee_arithmetic
    use fstats_distributions, only : distribution, probability_log, softplus, pi => distribution_pi
    use fstats_special_functions
    implicit none
    private
    public :: t_distribution

    type, extends(distribution) :: t_distribution
        !! Defines Student's T-Distribution.
        real(real64) :: dof
            !! The number of degrees of freedom.
    contains
        procedure, public :: log_pdf => td_log_pdf
        procedure, public :: survival => td_survival
        procedure, public :: log_cdf => td_log_cdf
        procedure, public :: log_survival => td_log_survival
        procedure, public :: pdf => td_pdf
        procedure, public :: cdf => td_cdf
        procedure, public :: mean => td_mean
        procedure, public :: median => td_median
        procedure, public :: mode => td_mode
        procedure, public :: variance => td_variance
        procedure, public :: recenter => td_recenter
    end type

contains

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
    rst = exp(this%log_pdf(x))
end function

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
    rst = this%survival(-x)
end function

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

pure function td_median(this) result(rst)
    !! Computes the median of the distribution.
    class(t_distribution), intent(in) :: this
        !! The t_distribution object.
    real(real64) :: rst

    ! Process
    rst = 0.0d0
end function

pure function td_mode(this) result(rst)
    !! Computes the mode of the distribution.
    class(t_distribution), intent(in) :: this
        !! The t_distribution object.
    real(real64) :: rst
        !! The mode.

    ! Process
    rst = 0.0d0
end function

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

subroutine td_recenter(this, x)
    !! Recenters the distribution about the supplied value.  This routine has
    !! no effect for this distribution as it is always centered about 0.
    class(t_distribution), intent(inout) :: this
        !! The t_distribution object.
    real(real64), intent(in) :: x
        !! The value about which to recenter.
end subroutine

pure elemental function td_log_pdf(this, x) result(rst)
    class(t_distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point; discrete masses require integer-valued x.
    real(real64) :: rst
        !! Log density; outside support gives -infinity, invalid inputs give NaN.
        !! An integrable density singularity at a support endpoint gives +infinity.
    real(real64) :: argument

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (ieee_is_nan(x)) return
    if (.not.ieee_is_finite(this%dof) .or. this%dof <= 0.0d0) return
    argument = 0.0d0
    if (x /= 0.0d0) argument = softplus(2.0d0 * log(abs(x)) - log(this%dof))
    rst = log_gamma(0.5d0 * this%dof + 0.5d0) - log_gamma(0.5d0 * this%dof) &
        - 0.5d0 * (log(this%dof) + log(pi)) - (0.5d0 * this%dof + 0.5d0) * argument
end function

pure elemental function td_survival(this, x) result(rst)
    class(t_distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point; discrete tails include masses strictly above x.
    real(real64) :: rst
        !! Probability on [0, 1]; invalid parameters or NaN x give NaN.
    real(real64) :: argument, tail

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (ieee_is_nan(x)) return
    if (.not.ieee_is_finite(this%dof) .or. this%dof <= 0.0d0) return
    argument = 1.0d0
    if (x /= 0.0d0) argument = exp(-softplus(2.0d0 * log(abs(x)) - log(this%dof)))
    tail = 0.5d0 * regularized_beta(0.5d0 * this%dof, 0.5d0, argument)
    rst = tail
    if (x < 0.0d0) rst = 1.0d0 - tail
end function

pure elemental function td_log_cdf(this, x) result(rst)
    class(t_distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point.
    real(real64) :: rst
        !! Log probability; zero probability gives -infinity, invalid inputs NaN.

    rst = probability_log(this%cdf(x))
    if (ieee_is_nan(rst)) return
    rst = this%log_survival(-x)
end function

pure elemental function td_log_survival(this, x) result(rst)
    class(t_distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point.
    real(real64) :: rst
        !! Log probability; zero probability gives -infinity, invalid inputs NaN.
    real(real64) :: argument, tail

    rst = probability_log(this%survival(x))
    if (ieee_is_nan(rst)) return
    argument = 1.0d0
    if (x /= 0.0d0) argument = exp(-softplus(2.0d0 * log(abs(x)) - log(this%dof)))
    tail = log_regularized_beta(0.5d0 * this%dof, 0.5d0, argument) - log(2.0d0)
    rst = tail
    if (x < 0.0d0) rst = stable_log1p(-exp(tail))
end function

end module fstats_t_distribution

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

module fstats_poisson_distribution
    use iso_fortran_env
    use ieee_arithmetic
    use fstats_distributions, only : distribution, probability_log
    use fstats_special_functions
    implicit none
    private
    public :: poisson_distribution

    type, extends(distribution) :: poisson_distribution
        !! Defines a Poisson distribution.
        real(real64), public :: occrence_rate
            !! The rate of occurrences.
    contains
        procedure, public :: log_pdf => pd_log_pdf
        procedure, public :: survival => pd_survival
        procedure, public :: log_cdf => pd_log_cdf
        procedure, public :: log_survival => pd_log_survival
        procedure, public :: pdf => pd_pdf
        procedure, public :: cdf => pd_cdf
        procedure, public :: mean => pd_mean
        procedure, public :: median => pd_median
        procedure, public :: mode => pd_mode
        procedure, public :: variance => pd_variance
        procedure, public :: recenter => pd_recenter
        procedure, public :: defined_range => pd_range
    end type

contains

pure elemental function pd_pdf(this, x) result(rst)
    !! Computes the probability mass function.
    !!
    !! The PMF for the Poisson distribution is given as \( f(k) = 
    !! \frac{\lambda^{k} e^{-\lambda}}{k!} \).
    class(poisson_distribution), intent(in) :: this
        !! The poisson_distribution object.
    real(real64), intent(in) :: x
        !! The number of occurrences (\(k\)).
    real(real64) :: rst
        !! The value of the function.
    
    ! Local Variables
    real(real64) :: lambda

    ! Process
    lambda = this%occrence_rate
    rst = exp(this%log_pdf(x))
end function

pure elemental function pd_cdf(this, x) result(rst)
    !! Computes the cumulative distribution function.
    !!
    !! The CDF for the Poisson distribution is given as 
    !! $$ F(k) = \frac{\Gamma(floor(k + 1), \lambda)}{floor(\lambda}!)} $$.
    class(poisson_distribution), intent(in) :: this
        !! The poisson_distribution object.
    real(real64), intent(in) :: x
        !! The number of occurrences (\(k\)).
    real(real64) :: rst
        !! The value of the function.

    ! Local Variables
    real(real64) :: lambda

    ! Process
    lambda = this%occrence_rate
    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (.not.ieee_is_finite(lambda) .or. ieee_is_nan(x)) return
    if (lambda < 0.0d0) return
    if (x < 0.0d0) then
        rst = 0.0d0
    else if (.not.ieee_is_finite(x)) then
        rst = 1.0d0
    else
        rst = regularized_gamma_upper(aint(x) + 1.0d0, lambda)
    end if
end function

pure function pd_mean(this) result(rst)
    !! Computes the mean of the distribution.
    class(poisson_distribution), intent(in) :: this
        !! The poisson_distribution object.
    real(real64) :: rst
        !! The mean.

    ! Process
    rst = this%occrence_rate
end function

pure function pd_median(this) result(rst)
    !! Computes the median of the distribution.
    class(poisson_distribution), intent(in) :: this
        !! The poisson_distribution object.
    real(real64) :: rst
        !! The median.

    real(real64) :: lambda
    lambda = this%occrence_rate
    rst = floor(lambda + 1.0d0 / 3.0d0 - 1.0d0 / (5.0d1 * lambda))
end function

pure function pd_mode(this) result(rst)
    !! Computes the mode of the distribution.
    class(poisson_distribution), intent(in) :: this
    !! The poisson_distribution object.
    real(real64) :: rst
        !! The mode.

    rst = max( &
        real(ceiling(this%occrence_rate) - 1.0d0, real64), &
        real(floor(this%occrence_rate), real64))
end function

pure function pd_variance(this) result(rst)
    !! Computes the variance.
    class(poisson_distribution), intent(in) :: this
        !! The poisson_distribution object.
    real(real64) :: rst
        !! The variance.

    rst = this%occrence_rate
end function

subroutine pd_recenter(this, x)
    !! Recenters the distribution about the supplied value.  This routine has
    !! no effect for this distribution as it is always centered about 0.
    class(poisson_distribution), intent(inout) :: this
        !! The poisson_distribution object.
    real(real64), intent(in) :: x
        !! The value about which to recenter.

    this%occrence_rate = x
end subroutine

pure function pd_range(this) result(rst)
    !! Gets the defined range for the distribution.
    class(poisson_distribution), intent(in) :: this
        !! The poisson_distribution object.
    real(real64), dimension(2) :: rst
        !! The defined range of the probability distributions [0, infinity).  As
        !! using a value of infinity may cause issue, this routine returns
        !! huge(0.0d0) instead.

    rst = [0.0d0, huge(0.0d0)]
end function

pure elemental function pd_log_pdf(this, x) result(rst)
    class(poisson_distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point; discrete masses require integer-valued x.
    real(real64) :: rst
        !! Log density; outside support gives -infinity, invalid inputs give NaN.
        !! An integrable density singularity at a support endpoint gives +infinity.

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (ieee_is_nan(x)) return
    if (.not.ieee_is_finite(this%occrence_rate) .or. this%occrence_rate < 0.0d0) return
    rst = ieee_value(0.0d0, ieee_negative_inf)
    if (x < 0.0d0 .or. .not.ieee_is_finite(x)) return
    if (x /= aint(x)) return
    if (this%occrence_rate == 0.0d0) then
        if (x == 0.0d0) rst = 0.0d0
    else
        rst = x * log(this%occrence_rate) - this%occrence_rate - log_gamma(x + 1.0d0)
    end if
end function

pure elemental function pd_survival(this, x) result(rst)
    class(poisson_distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point; discrete tails include masses strictly above x.
    real(real64) :: rst
        !! Probability on [0, 1]; invalid parameters or NaN x give NaN.

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (ieee_is_nan(x)) return
    if (.not.ieee_is_finite(this%occrence_rate) .or. this%occrence_rate < 0.0d0) return
    if (x < 0.0d0) then
        rst = 1.0d0
    else if (.not.ieee_is_finite(x)) then
        rst = 0.0d0
    else
        rst = regularized_gamma_lower(aint(x) + 1.0d0, this%occrence_rate)
    end if
end function

pure elemental function pd_log_cdf(this, x) result(rst)
    class(poisson_distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point.
    real(real64) :: rst
        !! Log probability; zero probability gives -infinity, invalid inputs NaN.

    rst = probability_log(this%cdf(x))
    if (ieee_is_nan(rst)) return
    if (x >= 0.0d0 .and. ieee_is_finite(x)) &
        rst = log_regularized_gamma_upper(aint(x) + 1.0d0, this%occrence_rate)
end function

pure elemental function pd_log_survival(this, x) result(rst)
    class(poisson_distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point.
    real(real64) :: rst
        !! Log probability; zero probability gives -infinity, invalid inputs NaN.

    rst = probability_log(this%survival(x))
    if (ieee_is_nan(rst)) return
    if (x >= 0.0d0 .and. ieee_is_finite(x)) &
        rst = log_regularized_gamma_lower(aint(x) + 1.0d0, this%occrence_rate)
end function

end module fstats_poisson_distribution

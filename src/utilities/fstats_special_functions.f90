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

module fstats_special_functions
    use iso_fortran_env
    use ieee_arithmetic
    implicit none
    private
    public :: beta
    public :: regularized_beta
    public :: incomplete_beta
    public :: incomplete_gamma_lower
    public :: incomplete_gamma_upper
    public :: regularized_gamma_lower
    public :: regularized_gamma_upper
    public :: digamma
    public :: log_regularized_beta
    public :: log_regularized_gamma_lower
    public :: log_regularized_gamma_upper
    public :: stable_log1p
    public :: stable_expm1

contains
pure elemental function stable_log1p(x) result(rst)
    !! Computes log(1 + x) accurately when x is near zero.
    real(real64), intent(in) :: x
        !! Argument on [-1, infinity); NaN inputs propagate.
    real(real64) :: rst
        !! Logarithm; x = -1 gives -infinity, x < -1 gives NaN.
    real(real64) :: shifted

    if (x == -1.0d0) then
        rst = ieee_value(0.0d0, ieee_negative_inf)
    else if (x < -1.0d0 .or. ieee_is_nan(x)) then
        rst = ieee_value(0.0d0, ieee_quiet_nan)
    else
        shifted = 1.0d0 + x
        if (shifted == 1.0d0) then
            rst = x
        else if (.not.ieee_is_finite(shifted)) then
            rst = log(shifted)
        else
            rst = log(shifted) * (x / (shifted - 1.0d0))
        end if
    end if
end function

pure elemental function stable_expm1(x) result(rst)
    !! Computes exp(x) - 1 without subtractive cancellation near zero.
    real(real64), intent(in) :: x
        !! Exponent argument.
    real(real64) :: rst
        !! Difference; overflow gives +infinity and NaN inputs propagate.
    real(real64) :: term
    integer(int32) :: iteration

    if (abs(x) >= 0.5d0 .or. ieee_is_nan(x)) then
        if (x > log(huge(1.0d0))) then
            rst = ieee_value(0.0d0, ieee_positive_inf)
        else
            rst = exp(x) - 1.0d0
        end if
    else
        term = x
        rst = x
        do iteration = 2, 100
            term = term * x / real(iteration, real64)
            rst = rst + term
            if (abs(term) <= epsilon(rst) * abs(rst)) exit
        end do
    end if
end function

pure elemental function log_complement(log_probability) result(rst)
    !! Computes log(1 - exp(log_probability)) without cancellation.
    real(real64), intent(in) :: log_probability
        !! Log probability on [-infinity, 0].
    real(real64) :: rst
        !! Log complement; zero input gives -infinity, positive/NaN inputs NaN.

    if (log_probability > 0.0d0 .or. ieee_is_nan(log_probability)) then
        rst = ieee_value(0.0d0, ieee_quiet_nan)
    else if (log_probability < -log(2.0d0)) then
        rst = stable_log1p(-exp(log_probability))
    else if (log_probability == 0.0d0) then
        rst = ieee_value(0.0d0, ieee_negative_inf)
    else
        rst = log(-stable_expm1(log_probability))
    end if
end function

pure elemental function log_regularized_beta(a, b, x, max_iterations) result(rst)
    !! Log regularized beta; invalid inputs or nonconvergence return NaN.
    !! The leading factor and complementary branch are evaluated in log space.
    real(real64), intent(in) :: a, b
        !! Finite, positive beta shape parameters.
    real(real64), intent(in) :: x
        !! Integration limit on [0, 1].
    integer(int32), intent(in), optional :: max_iterations
        !! Continued-fraction iteration limit; default 1000.
    real(real64) :: rst
        !! Log ratio; x = 0 gives -infinity, x = 1 gives zero.
    real(real64) :: leading, fraction

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (.not.ieee_is_finite(a) .or. .not.ieee_is_finite(b) .or. ieee_is_nan(x)) return
    if (a <= 0.0d0 .or. b <= 0.0d0 .or. x < 0.0d0 .or. x > 1.0d0) return
    if (x == 0.0d0) then
        rst = ieee_value(0.0d0, ieee_negative_inf)
        return
    else if (x == 1.0d0) then
        rst = 0.0d0
        return
    end if
    leading = log_gamma(a + b) - log_gamma(a) - log_gamma(b) &
        + a * log(x) + b * stable_log1p(-x)
    if (x < (a + 1.0d0) / (a + b + 2.0d0)) then
        fraction = beta_continued_fraction(a, b, x, max_iterations)
        if (.not.ieee_is_finite(fraction) .or. fraction <= 0.0d0) return
        rst = min(0.0d0, leading + log(fraction) - log(a))
    else
        fraction = beta_continued_fraction(b, a, 1.0d0 - x, max_iterations)
        if (.not.ieee_is_finite(fraction) .or. fraction <= 0.0d0) return
        rst = log_complement(min(0.0d0, leading + log(fraction) - log(b)))
    end if
end function

pure elemental function log_regularized_gamma_lower(a, x, max_iterations) result(rst)
    !! Log lower gamma ratio; invalid inputs or nonconvergence return NaN.
    real(real64), intent(in) :: a
        !! Finite, positive gamma shape parameter.
    real(real64), intent(in) :: x
        !! Nonnegative integration limit; +infinity is accepted.
    integer(int32), intent(in), optional :: max_iterations
        !! Series/continued-fraction iteration limit; default 10000.
    real(real64) :: rst
        !! Log P(a, x); x = 0 gives -infinity, x = +infinity gives zero.

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (.not.ieee_is_finite(a) .or. ieee_is_nan(x)) return
    if (a <= 0.0d0 .or. x < 0.0d0) return
    if (x == 0.0d0) then
        rst = ieee_value(0.0d0, ieee_negative_inf)
    else if (.not.ieee_is_finite(x)) then
        rst = 0.0d0
    else if (x < a + 1.0d0) then
        rst = log_gamma_series(a, x, max_iterations)
    else
        rst = log_complement(log_gamma_continued_fraction(a, x, max_iterations))
    end if
end function

pure elemental function log_regularized_gamma_upper(a, x, max_iterations) result(rst)
    !! Log upper gamma ratio; invalid inputs or nonconvergence return NaN.
    real(real64), intent(in) :: a
        !! Finite, positive gamma shape parameter.
    real(real64), intent(in) :: x
        !! Nonnegative integration limit; +infinity is accepted.
    integer(int32), intent(in), optional :: max_iterations
        !! Series/continued-fraction iteration limit; default 10000.
    real(real64) :: rst
        !! Log Q(a, x); x = 0 gives zero, x = +infinity gives -infinity.

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (.not.ieee_is_finite(a) .or. ieee_is_nan(x)) return
    if (a <= 0.0d0 .or. x < 0.0d0) return
    if (x == 0.0d0) then
        rst = 0.0d0
    else if (.not.ieee_is_finite(x)) then
        rst = ieee_value(0.0d0, ieee_negative_inf)
    else if (x < a + 1.0d0) then
        rst = log_complement(log_gamma_series(a, x, max_iterations))
    else
        rst = log_gamma_continued_fraction(a, x, max_iterations)
    end if
end function

! ------------------------------------------------------------------------------
pure elemental function beta(a, b) result(rst)
    !! Computes the beta function.
    !!
    !! The beta function is related to the gamma function
    !! by the following relationship.
    !! $$ \beta(a,b) = \frac{\Gamma(a) \Gamma(b)}{\Gamma(a + b)} $$.
    !!
    !! See Also
    !!
    !! - <a href="https://en.wikipedia.org/wiki/Beta_function" target="_blank">Wikipedia</a>
    real(real64), intent(in) :: a
        !! The first argument of the function.
    real(real64), intent(in) :: b
        !! The second argument of the function.
    real(real64) :: rst
        !! The value of the beta function at \( a \) and \( b \).

    ! Process
    ! REF: https://en.wikipedia.org/wiki/Beta_function
    ! LOG_GAMMA supplies the logarithm of the magnitude of the gamma function,
    ! so the sign must be restored separately for negative arguments.
    rst = gamma_sign(a) * gamma_sign(b) * gamma_sign(a + b) * &
        exp(log_gamma(a) + log_gamma(b) - log_gamma(a + b))
end function

! ------------------------------------------------------------------------------
pure elemental function gamma_sign(x) result(rst)
    !! Returns the sign of the gamma function.  The gamma function alternates
    !! in sign across successive intervals of unit width along the negative
    !! real axis.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the sign.
    real(real64) :: rst
        !! Either 1 or -1.

    ! Local Variables
    integer(int32) :: n

    if (x > 0.0d0) then
        rst = 1.0d0
    else
        n = int(floor(-x), int32)
        if (mod(n, 2) == 0) then
            rst = -1.0d0
        else
            rst = 1.0d0
        end if
    end if
end function

! ------------------------------------------------------------------------------
pure elemental function regularized_beta(a, b, x, max_iterations) result(rst)
    !! Computes the regularized beta function.
    !!
    !! The regularized beta function is defined as the ratio between
    !! the incomplete beta function and the beta function.
    !! $$ I_{x}(a,b) = \frac{\beta(x;a,b)}{\beta(a,b)} $$.
    !!
    !! Remarks
    !!
    !! The routine employs the continued fraction representation of the
    !! function, evaluated by means of the modified Lentz algorithm.  The
    !! leading factor is formed in logarithmic space such that the routine
    !! remains well-behaved for large arguments, where the beta function
    !! itself underflows.
    !!
    !! See Also
    !!
    !! - <a href="https://en.wikipedia.org/wiki/Beta_function" target="_blank">Wikipedia</a>
    real(real64), intent(in) :: a
        !! The first argument of the function.
    real(real64), intent(in) :: b
        !! The second argument of the function.
    real(real64), intent(in) :: x
        !! The upper limit of the integration.
    real(real64) :: rst
        !! The value of the regularized beta function.
    integer(int32), intent(in), optional :: max_iterations
        !! Iteration limit; nonconvergence returns NaN.

    rst = exp(log_regularized_beta(a, b, x, max_iterations))
end function

! ------------------------------------------------------------------------------
pure elemental function beta_continued_fraction(a, b, x, max_iterations) result(rst)
    !! Evaluates the continued fraction expansion of the incomplete beta
    !! function by means of the modified Lentz algorithm.
    real(real64), intent(in) :: a
        !! The first argument of the function.
    real(real64), intent(in) :: b
        !! The second argument of the function.
    real(real64), intent(in) :: x
        !! The upper limit of the integration.
    real(real64) :: rst
        !! The value of the continued fraction.
    integer(int32), intent(in), optional :: max_iterations
        !! Iteration limit, default 1000; exhaustion returns NaN.

    ! Parameters
    integer(int32), parameter :: maxiter = 1000
    real(real64), parameter :: one = 1.0d0
    real(real64), parameter :: two = 2.0d0

    ! Local Variables
    integer(int32) :: m, m2, limit
    real(real64) :: aa, c, d, del, qab, qam, qap, eps, fpmin, rm

    ! Initialization
    eps = epsilon(eps)
    limit = maxiter
    if (present(max_iterations)) limit = max_iterations
    fpmin = tiny(fpmin) / eps
    qab = a + b
    qap = a + one
    qam = a - one
    c = one
    d = one - qab * x / qap
    if (abs(d) < fpmin) d = sign(fpmin, d)
    d = one / d
    rst = d

    ! Process
    do m = 1, limit
        rm = real(m, real64)
        m2 = 2 * m

        ! The even step of the recurrence
        aa = rm * (b - rm) * x / ((qam + m2) * (a + m2))
        d = one + aa * d
        if (abs(d) < fpmin) d = sign(fpmin, d)
        c = one + aa / c
        if (abs(c) < fpmin) c = sign(fpmin, c)
        d = one / d
        rst = rst * d * c

        ! The odd step of the recurrence
        aa = -(a + rm) * (qab + rm) * x / ((a + m2) * (qap + m2))
        d = one + aa * d
        if (abs(d) < fpmin) d = sign(fpmin, d)
        c = one + aa / c
        if (abs(c) < fpmin) c = sign(fpmin, c)
        d = one / d
        del = d * c
        rst = rst * del

        if (abs(del - one) <= eps) return
    end do
    rst = ieee_value(0.0d0, ieee_quiet_nan)
end function

! ------------------------------------------------------------------------------
pure elemental function incomplete_beta(a, b, x) result(rst)
    !! Computes the incomplete beta function.
    !!
    !! The incomplete beta function is defind as:
    !! $$ \beta(x;a,b) = \int_{0}^{x} t^{a-1} (1 - t)^{b-1} dt $$.
    !!
    !! See Also
    !!
    !! - <a href="https://en.wikipedia.org/wiki/Beta_function#Incomplete_beta_function" target="_blank">Wikipedia</a>
    real(real64), intent(in) :: a
        !! The first argument of the function.
    real(real64), intent(in) :: b
        !! The second argument of the function.
    real(real64), intent(in) :: x
        !! The upper limit of the integration.
    real(real64) :: rst
        !! The value of the incomplete beta function.

    ! Process
    rst = beta(a, b) * regularized_beta(a, b, x)
end function

! ------------------------------------------------------------------------------
pure elemental function regularized_gamma_lower(a, x, max_iterations) result(rst)
    !! Computes the regularized lower incomplete gamma function.
    !!
    !! The regularized lower incomplete gamma function is defined as:
    !! $$ P(a, x) = \frac{\gamma(a, x)}{\Gamma(a)} $$
    !!
    !! Remarks
    !!
    !! The function is evaluated by means of a series expansion for 
    !! \( x < a + 1 \) and a continued fraction otherwise, both carried to
    !! convergence and formed in logarithmic space.  The routine is therefore
    !! well-behaved for arbitrarily large arguments, unlike the unregularized
    !! forms whose values overflow.
    !!
    !! See Also
    !!
    !! - <a href="https://en.wikipedia.org/wiki/Incomplete_gamma_function" target="_blank">Wikipedia</a>
    real(real64), intent(in) :: a
        !! The coefficient value.  The value must be positive.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.  The value must be
        !! non-negative.
    real(real64) :: rst
        !! The function value, which lies on the interval [0, 1].
    integer(int32), intent(in), optional :: max_iterations
        !! Iteration limit; nonconvergence returns NaN.

    rst = exp(log_regularized_gamma_lower(a, x, max_iterations))
end function

! ------------------------------------------------------------------------------
pure elemental function regularized_gamma_upper(a, x, max_iterations) result(rst)
    !! Computes the regularized upper incomplete gamma function.
    !!
    !! The regularized upper incomplete gamma function is defined as:
    !! $$ Q(a, x) = \frac{\Gamma(a, x)}{\Gamma(a)} = 1 - P(a, x) $$
    !!
    !! Remarks
    !!
    !! Whichever of the series or the continued fraction converges rapidly is
    !! evaluated directly, such that the result retains its relative accuracy
    !! even when the probability is vanishingly small.
    !!
    !! See Also
    !!
    !! - <a href="https://en.wikipedia.org/wiki/Incomplete_gamma_function" target="_blank">Wikipedia</a>
    real(real64), intent(in) :: a
        !! The coefficient value.  The value must be positive.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.  The value must be
        !! non-negative.
    real(real64) :: rst
        !! The function value, which lies on the interval [0, 1].
    integer(int32), intent(in), optional :: max_iterations
        !! Iteration limit; nonconvergence returns NaN.

    rst = exp(log_regularized_gamma_upper(a, x, max_iterations))
end function

! ------------------------------------------------------------------------------
pure elemental function log_gamma_series(a, x, max_iterations) result(rst)
    !! Evaluates the log lower gamma ratio by a convergent series for x < a + 1.
    real(real64), intent(in) :: a
        !! The coefficient value.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! Log ratio; nonconvergence returns NaN.
    integer(int32), intent(in), optional :: max_iterations
        !! Iteration limit, default 10000.

    ! Parameters
    integer(int32), parameter :: maxiter = 10000

    ! Local Variables
    integer(int32) :: i, limit
    real(real64) :: ap, del, s, eps

    eps = epsilon(eps)
    limit = maxiter
    if (present(max_iterations)) limit = max_iterations
    ap = a
    s = 1.0d0 / a
    del = s
    rst = ieee_value(0.0d0, ieee_quiet_nan)
    do i = 1, limit
        ap = ap + 1.0d0
        del = del * x / ap
        s = s + del
        if (abs(del) < abs(s) * eps) then
            rst = min(0.0d0, log(s) - x + a * log(x) - log_gamma(a))
            return
        end if
    end do
end function

! ------------------------------------------------------------------------------
pure elemental function log_gamma_continued_fraction(a, x, max_iterations) result(rst)
    !! Evaluates the log upper gamma ratio using modified Lentz iteration.
    !! The expansion converges rapidly for x >= a + 1.
    real(real64), intent(in) :: a
        !! The coefficient value.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! Log ratio; nonconvergence returns NaN.
    integer(int32), intent(in), optional :: max_iterations
        !! Iteration limit, default 10000.

    ! Parameters
    integer(int32), parameter :: maxiter = 10000

    ! Local Variables
    integer(int32) :: i, limit
    real(real64) :: an, b, c, d, del, h, eps, fpmin

    eps = epsilon(eps)
    limit = maxiter
    if (present(max_iterations)) limit = max_iterations
    fpmin = tiny(fpmin) / eps
    b = x + 1.0d0 - a
    c = 1.0d0 / fpmin
    d = 1.0d0 / b
    h = d
    rst = ieee_value(0.0d0, ieee_quiet_nan)
    do i = 1, limit
        an = -i * (i - a)
        b = b + 2.0d0
        d = an * d + b
        if (abs(d) < fpmin) d = sign(fpmin, d)
        c = b + an / c
        if (abs(c) < fpmin) c = sign(fpmin, c)
        d = 1.0d0 / d
        del = d * c
        h = h * del
        if (abs(del - 1.0d0) <= eps) then
            if (h > 0.0d0 .and. ieee_is_finite(h)) &
                rst = min(0.0d0, log(h) - x + a * log(x) - log_gamma(a))
            return
        end if
    end do
end function

! ------------------------------------------------------------------------------
pure elemental function incomplete_gamma_upper(a, x) result(rst)
    !! Computes the upper incomplete gamma function.
    !!
    !! The upper incomplete gamma function is defined as:
    !! $$ \Gamma(a, x) = \int_{x}^{\infty} t^{a-1} e^{-t} \,dt $$
    !!
    !! Remarks
    !!
    !! The value overflows for arguments beyond roughly \( a = 171 \), where
    !! \( \Gamma(a) \) itself exceeds the range of a double precision number.
    !! Use [[regularized_gamma_upper]] where only the ratio to \( \Gamma(a) \)
    !! is required.
    !!
    !! See Also
    !!
    !! - <a href="https://en.wikipedia.org/wiki/Incomplete_gamma_function" target="_blank">Wikipedia</a>
    real(real64), intent(in) :: a
        !! The coefficient value.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The function value.

    rst = regularized_gamma_upper(a, x) * gamma(a)
end function

! ------------------------------------------------------------------------------
pure elemental function incomplete_gamma_lower(a, x) result(rst)
    !! Computes the lower incomplete gamma function.
    !!
    !! The lower incomplete gamma function is defined as:
    !! $$ \gamma(a, x) = \int_{0}^{x} t^{a-1} e^{-t} \,dt $$
    !!
    !! Remarks
    !!
    !! The value overflows for arguments beyond roughly \( a = 171 \), where
    !! \( \Gamma(a) \) itself exceeds the range of a double precision number.
    !! Use [[regularized_gamma_lower]] where only the ratio to \( \Gamma(a) \)
    !! is required.
    !!
    !! See Also
    !!
    !! - <a href="https://en.wikipedia.org/wiki/Incomplete_gamma_function" target="_blank">Wikipedia</a>
    real(real64), intent(in) :: a
        !! The coefficient value.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The function value.

    rst = regularized_gamma_lower(a, x) * gamma(a)
end function

! ------------------------------------------------------------------------------
pure elemental function digamma(x) result(rst)
    !! Computes the digamma function.
    !!
    !! The digamma function is defined as:
    !! $$ \psi(x) = 
    !! \frac{d}{dx}\left( \ln \left( \Gamma \left( x \right) \right) 
    !! \right) $$
    !!
    !! See Also
    !!
    !! - <a href="https://en.wikipedia.org/wiki/Digamma_function" target="_blank">Wikipedia</a>
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The function value.

    ! Parameters
    ! The asymptotic expansion below is truncated after the 1/x**10 term, so
    ! the threshold must be large enough that the first neglected term,
    ! 691 / (32760 x**12), falls below the rounding level.
    real(real64), parameter :: c = 2.0d1
    real(real64), parameter :: euler_mascheroni = 0.57721566490153286060d0

    ! Local Variables
    real(real64) :: r, x2, nan
    
    ! REF:
    ! - https://people.sc.fsu.edu/~jburkardt/f_src/asa103/asa103.f90

    ! If x <= 0.0
    if (x <= 0.0) then
        nan = ieee_value(nan, IEEE_QUIET_NAN)
        rst = nan
        return
    end if

    ! Approximation for a small argument
    if (x <= 1.0d-6) then
        rst = -euler_mascheroni - 1.0d0 / x + 1.6449340668482264365d0 * x
        return
    end if

    ! Process
    rst = 0.0d0
    x2 = x
    do while (x2 < c)
        rst = rst - 1.0d0 / x2
        x2 = x2 + 1.0d0
    end do

    r = 1.0d0 / x2
    rst = rst + log(x2) - 0.5d0 * r
    r = r * r
    rst = rst &
        -r * (1.0d0 / 12.0d0 &
        - r * (1.0d0 / 120.0d0 &
        - r * (1.0d0 / 252.0d0 &
        - r * (1.0d0 / 240.0d0 &
        - r * (1.0d0 / 132.0d0) &
    ))))
end function

! ------------------------------------------------------------------------------
end module
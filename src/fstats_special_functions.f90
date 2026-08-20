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

contains
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
pure elemental function regularized_beta(a, b, x) result(rst)
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

    ! Local Variables
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0
    real(real64), parameter :: two = 2.0d0
    real(real64) :: bt

    ! Handle the limits of the interval directly
    if (x <= zero) then
        rst = zero
        return
    end if
    if (x >= one) then
        rst = one
        return
    end if

    ! The leading factor x**a * (1 - x)**b / beta(a,b)
    bt = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + &
        a * log(x) + b * log(one - x))

    ! The continued fraction converges rapidly only for x below the mean of
    ! the distribution; the symmetry relation handles the remainder
    if (x < (a + one) / (a + b + two)) then
        rst = bt * beta_continued_fraction(a, b, x) / a
    else
        rst = one - bt * beta_continued_fraction(b, a, one - x) / b
    end if
end function

! ------------------------------------------------------------------------------
pure elemental function beta_continued_fraction(a, b, x) result(rst)
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

    ! Parameters
    integer(int32), parameter :: maxiter = 1000
    real(real64), parameter :: one = 1.0d0
    real(real64), parameter :: two = 2.0d0

    ! Local Variables
    integer(int32) :: m, m2
    real(real64) :: aa, c, d, del, qab, qam, qap, eps, fpmin, rm

    ! Initialization
    eps = epsilon(eps)
    fpmin = tiny(fpmin) / eps
    qab = a + b
    qap = a + one
    qam = a - one
    c = one
    d = one - qab * x / qap
    if (abs(d) < fpmin) d = fpmin
    d = one / d
    rst = d

    ! Process
    do m = 1, maxiter
        rm = real(m, real64)
        m2 = 2 * m

        ! The even step of the recurrence
        aa = rm * (b - rm) * x / ((qam + m2) * (a + m2))
        d = one + aa * d
        if (abs(d) < fpmin) d = fpmin
        c = one + aa / c
        if (abs(c) < fpmin) c = fpmin
        d = one / d
        rst = rst * d * c

        ! The odd step of the recurrence
        aa = -(a + rm) * (qab + rm) * x / ((a + m2) * (qap + m2))
        d = one + aa * d
        if (abs(d) < fpmin) d = fpmin
        c = one + aa / c
        if (abs(c) < fpmin) c = fpmin
        d = one / d
        del = d * c
        rst = rst * del

        if (abs(del - one) <= eps) exit
    end do
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
pure elemental function regularized_gamma_lower(a, x) result(rst)
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

    if (a <= 0.0d0 .or. x < 0.0d0) then
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    else if (x == 0.0d0) then
        rst = 0.0d0
    else if (x < a + 1.0d0) then
        rst = gamma_series(a, x)
    else
        rst = 1.0d0 - gamma_continued_fraction(a, x)
    end if
end function

! ------------------------------------------------------------------------------
pure elemental function regularized_gamma_upper(a, x) result(rst)
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

    if (a <= 0.0d0 .or. x < 0.0d0) then
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    else if (x == 0.0d0) then
        rst = 1.0d0
    else if (x < a + 1.0d0) then
        rst = 1.0d0 - gamma_series(a, x)
    else
        rst = gamma_continued_fraction(a, x)
    end if
end function

! ------------------------------------------------------------------------------
pure elemental function gamma_series(a, x) result(rst)
    !! Evaluates the series representation of the regularized lower incomplete
    !! gamma function.  The series converges rapidly for x below a + 1.
    real(real64), intent(in) :: a
        !! The coefficient value.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The function value.

    ! Parameters
    integer(int32), parameter :: maxiter = 10000

    ! Local Variables
    integer(int32) :: i
    real(real64) :: ap, del, s, eps

    eps = epsilon(eps)
    ap = a
    s = 1.0d0 / a
    del = s
    do i = 1, maxiter
        ap = ap + 1.0d0
        del = del * x / ap
        s = s + del
        if (abs(del) < abs(s) * eps) exit
    end do
    rst = s * exp(-x + a * log(x) - log_gamma(a))
end function

! ------------------------------------------------------------------------------
pure elemental function gamma_continued_fraction(a, x) result(rst)
    !! Evaluates the continued fraction representation of the regularized 
    !! upper incomplete gamma function by means of the modified Lentz 
    !! algorithm.  The expansion converges rapidly for x above a + 1.
    real(real64), intent(in) :: a
        !! The coefficient value.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The function value.

    ! Parameters
    integer(int32), parameter :: maxiter = 10000

    ! Local Variables
    integer(int32) :: i
    real(real64) :: an, b, c, d, del, h, eps, fpmin

    eps = epsilon(eps)
    fpmin = tiny(fpmin) / eps
    b = x + 1.0d0 - a
    c = 1.0d0 / fpmin
    d = 1.0d0 / b
    h = d
    do i = 1, maxiter
        an = -i * (i - a)
        b = b + 2.0d0
        d = an * d + b
        if (abs(d) < fpmin) d = fpmin
        c = b + an / c
        if (abs(c) < fpmin) c = fpmin
        d = 1.0d0 / d
        del = d * c
        h = h * del
        if (abs(del - 1.0d0) <= eps) exit
    end do
    rst = h * exp(-x + a * log(x) - log_gamma(a))
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
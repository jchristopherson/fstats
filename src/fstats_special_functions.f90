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
    rst = exp(log_gamma(a) + log_gamma(b) - log_gamma(a + b))
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
! REF: https://people.math.sc.edu/Burkardt/f_src/special_functions/special_functions.f90
pure elemental function incomplete_gamma_upper(a, x) result(rst)
    !! Computes the upper incomplete gamma function.
    !!
    !! The upper incomplete gamma function is defined as:
    !! $$ \Gamma(a, x) = \int_{x}^{\infty} t^{a-1} e^{-t} \,dt $$
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

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0
    real(real64), parameter :: ten = 1.0d1

    ! Local Variables
    real(real64) :: ga, gin, gip, r, s, t0, xam, small
    integer(int32) :: k

    ! Process
    small = ten * epsilon(small)
    xam = -x + a * log(x)
    if (xam > 7.0d2 .or. a > 1.7d2) then
        rst = ieee_value(rst, IEEE_QUIET_NAN)
        return
    end if

    if (x == zero) then
        rst = gamma(a)
    else if (x <= one + a) then
        s = one / a
        r = s
        do k = 1, 60
            r = r * x / (a + k)
            s = s + r
            if (abs(r / s) < small) then
                exit
            end if
        end do

        gin = exp(xam) * s
        ga = gamma(a)
        gip = gin / ga
        rst = ga - gin
    else if (one + a < x) then
        t0 = zero
        do k = 60, 1, -1
            t0 = (k - a) / (one + k / (x + t0))
        end do
        rst = exp(xam) / (x + t0)
    end if
end function

! ------------------------------------------------------------------------------
pure elemental function incomplete_gamma_lower(a, x) result(rst)
    !! Computes the lower incomplete gamma function.
    !!
    !! The lower incomplete gamma function is defined as:
    !! $$ \gamma(a, x) = \int_{0}^{x} t^{a-1} e^{-t} \,dt $$
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

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0
    real(real64), parameter :: ten = 1.0d1

    ! Local Variables
    real(real64) :: ga, gim, r, s, t0, xam, small
    integer(int32) :: k

    ! Process
    small = ten * epsilon(small)
    xam = -x + a * log(x)
    if (xam > 7.0d2 .or. a > 1.7d2) then
        rst = ieee_value(rst, IEEE_QUIET_NAN)
        return
    end if

    if (x == zero) then
        rst = 0.0d0
    else if (x <= one + a) then
        s = one / a
        r = s
        do k = 1, 60
            r = r * x / (a + k)
            s = s + r
            if (abs(r / s) < small) then
                exit
            end if
        end do

        rst = exp(xam) * s
    else if (one + a < x) then
        t0 = zero
        do k = 60, 1, -1
            t0 = (k - a) / (one + k / (x + t0))
        end do
        gim = exp(xam) / (x + t0)
        ga = gamma(a)
        rst = ga - gim
    end if
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
    real(real64), parameter :: c = 8.5d0
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
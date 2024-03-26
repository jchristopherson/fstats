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
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Beta_function)
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
! source: https://people.math.sc.edu/Burkardt/f_src/special_functions/special_functions.f90
pure elemental function regularized_beta(a, b, x) result(rst)
    !! Computes the regularized beta function.
    !!
    !! The regularized beta function is defined as the ratio between
    !! the incomplete beta function and the beta function.
    !! $$ I_{x}(a,b) = \frac{\beta(x;a,b)}{\beta(a,b)} $$.
    !!
    !! See Also
    !!
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Beta_function)
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
    real(real64) :: bt, dk(51), fk(51), s0, t1, t2, ta, tb
    integer(int32) :: k

    ! Process
    s0 = (a + one) / (a + b + two)
    bt = beta(a, b)

    if (x <= s0) then
        do k = 1, 20
            dk(2*k) = k * (b - k) * x / (a + two * k - one) / (a + two * k)
        end do
        do k = 0, 20
            dk(2*k+1) = -(a + k) * (a + b + k) * x / (a + two * k) / &
                (a + two * k + one)
        end do

        t1 = zero
        do k = 20, 1, -1
            t1 = dk(k) / (one + t1)
        end do
        ta = one / (one + t1)
        rst = x**a * (one - x)**b / (a * bt) * ta
    else
        do k = 1, 20
            fk(2*k) = k * (a - k) * (one - x) / (b + two * k - one) / &
                (b + two * k)
        end do
        do k = 0, 20
            fk(2*k+1) = -(b + k) * (a + b + k) * (one - x) / (b + two * k) / &
                (b + two * k + one)
        end do

        t2 = zero
        do k = 20, 1, -1
            t2 = fk(k) / (one + t2)
        end do
        tb = one / (one + t2)
        rst = one - x**a * (one - x)**b / (b * bt) * tb
    end if
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
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Beta_function#Incomplete_beta_function)
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
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Incomplete_gamma_function)
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
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Incomplete_gamma_function)
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
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Digamma_function)
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
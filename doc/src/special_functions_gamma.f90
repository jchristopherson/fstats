submodule (fstats) special_functions_gamma
    use ieee_arithmetic
contains
! ------------------------------------------------------------------------------
! REF: https://people.math.sc.edu/Burkardt/f_src/special_functions/special_functions.f90
pure elemental module function incomplete_gamma_upper_real64(a, x) result(rst)
    ! Arguments
    real(real64), intent(in) :: a, x
    real(real64) :: rst

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
pure elemental module function incomplete_gamma_upper_real32(a, x) result(rst)
    ! Arguments
    real(real32), intent(in) :: a, x
    real(real32) :: rst

    ! Parameters
    real(real32), parameter :: zero = 0.0
    real(real32), parameter :: one = 1.0
    real(real32), parameter :: ten = 1.0e1

    ! Local Variables
    real(real32) :: ga, gin, gip, r, s, t0, xam, small
    integer(int32) :: k

    ! Process
    small = ten * epsilon(small)
    xam = -x + a * log(x)
    if (xam > 7.0e2 .or. a > 1.7e2) then
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
pure elemental module function incomplete_gamma_lower_real64(a, x) result(rst)
    ! Arguments
    real(real64), intent(in) :: a, x
    real(real64) :: rst

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
pure elemental module function incomplete_gamma_lower_real32(a, x) result(rst)
    ! Arguments
    real(real32), intent(in) :: a, x
    real(real32) :: rst

    ! Parameters
    real(real32), parameter :: zero = 0.0
    real(real32), parameter :: one = 1.0
    real(real32), parameter :: ten = 1.0e1

    ! Local Variables
    real(real32) :: ga, gim, r, s, t0, xam, small
    integer(int32) :: k

    ! Process
    small = ten * epsilon(small)
    xam = -x + a * log(x)
    if (xam > 7.0e2 .or. a > 1.7e2) then
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
end submodule
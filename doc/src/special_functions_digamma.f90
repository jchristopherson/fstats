submodule (fstats) special_functions_digamma
    use ieee_arithmetic
contains
! ------------------------------------------------------------------------------
pure elemental module function digamma_real64(x) result(rst)
    ! Arguments
    real(real64), intent(in) :: x
    real(real64) :: rst

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
pure elemental module function digamma_real32(x) result(rst)
    ! Arguments
    real(real32), intent(in) :: x
    real(real32) :: rst

    ! Parameters
    real(real32), parameter :: c = 8.5
    real(real32), parameter :: euler_mascheroni = 0.57721566490153286060

    ! Local Variables
    real(real32) :: r, x2, nan
    
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
        rst = -euler_mascheroni - 1.0 / x + 1.6449340668482264365 * x
        return
    end if

    ! Process
    rst = 0.0
    x2 = x
    do while (x2 < c)
        rst = rst - 1.0 / x2
        x2 = x2 + 1.0
    end do

    r = 1.0 / x2
    rst = rst + log(x2) - 0.5 * r
    r = r * r
    rst = rst &
        -r * (1.0 / 12.0 &
        - r * (1.0 / 120.0 &
        - r * (1.0 / 252.0 &
        - r * (1.0 / 240.0 &
        - r * (1.0 / 132.0) &
    ))))
end function

! ------------------------------------------------------------------------------
end submodule
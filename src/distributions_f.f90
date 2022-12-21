submodule (fstats) distributions_f
    use ieee_arithmetic
contains
! ------------------------------------------------------------------------------
pure module elemental function fd_pdf(this, x) result(rst)
    ! Arguments
    class(f_distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst

    ! Process
    real(real64) :: d1, d2
    d1 = this%d1
    d2 = this%d2
    rst = (1.0d0 / beta(0.5d0 * d1, 0.5d0 * d2)) * (d1 / d2)**(0.5d0 * d1) * &
        x**(0.5d0 * d1 - 1.0d0) * (1.0d0 + d1 * x/ d2)**(-0.5d0 * (d1 + d2))
end function

! ------------------------------------------------------------------------------
pure module elemental function fd_cdf(this, x) result(rst)
    ! Arguments
    class(f_distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst

    ! Process
    real(real64) :: d1, d2
    d1 = this%d1
    d2 = this%d2
    rst = regularized_beta(0.5d0 * d1, 0.5d0 * d2, d1 * x / (d1 * x + d2))
end function

! ------------------------------------------------------------------------------
pure module function fd_mean(this) result(rst)
    ! Arguments
    class(f_distribution), intent(in) :: this
    real(real64) :: rst

    ! Process
    if (this%d2 > 2.0d0) then
        rst = this%d2 / (this%d2 - 2.0d0)
    else
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    end if
end function

! ------------------------------------------------------------------------------
pure module function fd_median(this) result(rst)
    ! Arguments
    class(f_distribution), intent(in) :: this
    real(real64) :: rst
    rst = ieee_value(rst, IEEE_QUIET_NAN)
end function

! ------------------------------------------------------------------------------
pure module function fd_mode(this) result(rst)
    ! Arguments
    class(f_distribution), intent(in) :: this
    real(real64) :: rst

    ! Process
    if (this%d1 > 2.0d0) then
        rst = ((this%d1 - 2.0d0) / this%d1) * (this%d2 / (this%d2 + 2.0d0))
    else
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    end if
end function

! ------------------------------------------------------------------------------
pure module function fd_variance(this) result(rst)
    ! Arguments
    class(f_distribution), intent(in) :: this
    real(real64) :: rst

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

! ------------------------------------------------------------------------------
end submodule
submodule (fstats) distributions_t
    use ieee_arithmetic
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
contains
! ------------------------------------------------------------------------------
! REF: https://en.wikipedia.org/wiki/Student%27s_t-distribution
pure module elemental function td_pdf(this, x) result(rst)
    ! Arguments
    class(t_distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst

    ! Process
    rst = gamma((this%dof + 1.0d0) / 2.0d0) / &
        (sqrt(this%dof * pi) * gamma(this%dof / 2.0d0)) *&
        (1.0d0 + x**2 / this%dof)**(-0.5d0 * (1.0d0 + this%dof))
end function

! ------------------------------------------------------------------------------
pure module elemental function td_cdf(this, x) result(rst)
    ! Arguments
    class(t_distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst

    ! Process
    real(real64) :: t
    t = this%dof / (this%dof + x**2)
    rst = 1.0d0 - 0.5d0 * regularized_beta(0.5d0 * this%dof, 0.5d0, t)
    if (x < 0) rst = 1.0d0 - rst
end function

! ------------------------------------------------------------------------------
pure module function td_mean(this) result(rst)
    class(t_distribution), intent(in) :: this
    real(real64) :: rst
    if (this%dof < 1.0d0) then
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    else
        rst = 0.0d0
    end if
end function

! ------------------------------------------------------------------------------
pure module function td_median(this) result(rst)
    class(t_distribution), intent(in) :: this
    real(real64) :: rst
    rst = 0.0d0
end function

! ------------------------------------------------------------------------------
pure module function td_mode(this) result(rst)
    class(t_distribution), intent(in) :: this
    real(real64) :: rst
    rst = 0.0d0
end function

! ------------------------------------------------------------------------------
pure module function td_variance(this) result(rst)
    class(t_distribution), intent(in) :: this
    real(real64) :: rst
    if (this%dof <= 1.0d0) then
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    else if (this%dof > 1.0d0 .and. this%dof <= 2.0d0) then
        rst = ieee_value(rst, IEEE_POSITIVE_INF)
    else
        rst = this%dof / (this%dof - 2.0d0)
    end if
end function

! ------------------------------------------------------------------------------
end submodule
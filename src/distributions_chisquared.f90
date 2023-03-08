submodule (fstats) distributions_chisquared
contains
! ------------------------------------------------------------------------------
pure module elemental function cs_pdf(this, x) result(rst)
    ! Arguments
    class(chi_squared_distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst

    ! Local Variables
    real(real64) :: arg

    ! Process
    arg = 0.5d0 * this%dof
    rst = 1.0d0 / (2.0d0**arg * gamma(arg)) * x**(arg - 1.0d0) * exp(-0.5d0 * x)
end function

! ------------------------------------------------------------------------------
pure module elemental function cs_cdf(this, x) result(rst)
    ! Arguments
    class(chi_squared_distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst

    ! Local Variables
    real(real64) :: arg

    ! Process
    arg = 0.5d0 * this%dof
    rst = incomplete_gamma_lower(arg, 0.5d0 * x) / gamma(arg)
end function

! ------------------------------------------------------------------------------
pure module function cs_mean(this) result(rst)
    ! Arguments
    class(chi_squared_distribution), intent(in) :: this
    real(real64) :: rst

    ! Process
    rst = real(this%dof, real64)
end function

! ------------------------------------------------------------------------------
pure module function cs_median(this) result(rst)
    ! Arguments
    class(chi_squared_distribution), intent(in) :: this
    real(real64) :: rst

    ! Process
    rst = this%dof * (1.0d0 - 2.0d0 / (9.0d0 * this%dof))**3
end function

! ------------------------------------------------------------------------------
pure module function cs_mode(this) result(rst)
    ! Arguments
    class(chi_squared_distribution), intent(in) :: this
    real(real64) :: rst

    ! Process
    rst = max(this%dof - 2.0d0, 0.0d0)
end function

! ------------------------------------------------------------------------------
pure module function cs_variance(this) result(rst)
    ! Arguments
    class(chi_squared_distribution), intent(in) :: this
    real(real64) :: rst

    ! Process
    rst = 2.0d0 * this%dof
end function

! ------------------------------------------------------------------------------
end submodule
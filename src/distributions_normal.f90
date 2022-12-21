submodule (fstats) distributions_normal
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
contains
! ------------------------------------------------------------------------------
pure module elemental function nd_pdf(this, x) result(rst)
    class(normal_distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst
    rst = exp(-0.5d0 * ((x - this%mean_value) / this%standard_deviation)**2) / &
        (this%standard_deviation * sqrt(2.0d0 * pi))
end function

! ------------------------------------------------------------------------------
pure module elemental function nd_cdf(this, x) result(rst)
    class(normal_distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst
    rst = 0.5d0 * (1.0d0 + erf((x - this%mean_value) / &
        (this%standard_deviation * sqrt(2.0d0))))
end function

! ------------------------------------------------------------------------------
pure module function nd_mean(this) result(rst)
    class(normal_distribution), intent(in) :: this
    real(real64) :: rst
    rst = this%mean_value
end function

! ------------------------------------------------------------------------------
pure module function nd_median(this) result(rst)
    class(normal_distribution), intent(in) :: this
    real(real64) :: rst
    rst = this%mean_value
end function

! ------------------------------------------------------------------------------
pure module function nd_mode(this) result(rst)
    class(normal_distribution), intent(in) :: this
    real(real64) :: rst
    rst = this%mean_value
end function

! ------------------------------------------------------------------------------
pure module function nd_variance(this) result(rst)
    class(normal_distribution), intent(in) :: this
    real(real64) :: rst
    rst = this%standard_deviation**2
end function

! ------------------------------------------------------------------------------
module subroutine nd_standardize(this)
    class(normal_distribution), intent(inout) :: this
    this%mean_value = 0.0d0
    this%standard_deviation = 1.0d0
end subroutine

! ------------------------------------------------------------------------------
end submodule
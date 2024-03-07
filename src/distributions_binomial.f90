submodule (fstats) distributions_binomial
    implicit none
contains
! ------------------------------------------------------------------------------
pure elemental function fact(x) result(rst)
    !! Computes the factorial of X.
    real(real64), intent(in) :: x
    real(real64) :: rst
    rst = gamma(x + 1.0d0)
end function

! ------------------------------------------------------------------------------
pure module elemental function bd_pdf(this, x) result(rst)
    ! Arguments
    class(binomial_distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst

    ! Local Variables
    real(real64) :: dn

    ! Process
    dn = real(this%n, real64)
    rst = (fact(dn) / (fact(x) * fact(dn - x))) * (this%p**x) * (1.0d0 - this%p)**(dn - x)
end function

! ------------------------------------------------------------------------------
pure module elemental function bd_cdf(this, x) result(rst)
    ! Arguments
    class(binomial_distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst

    ! Local Variables
    real(real64) :: dn

    ! Process
    dn = real(this%n, real64)
    rst = regularized_beta(dn - x, x + 1.0d0, 1.0d0 - this%p)
end function

! ------------------------------------------------------------------------------
pure module function bd_mean(this) result(rst)
    class(binomial_distribution), intent(in) :: this
    real(real64) :: rst
    rst = real(this%n * this%p, real64)
end function

! ------------------------------------------------------------------------------
pure module function bd_median(this) result(rst)
    class(binomial_distribution), intent(in) :: this
    real(real64) :: rst
    rst = real(this%n * this%p, real64)
end function

! ------------------------------------------------------------------------------
pure module function bd_mode(this) result(rst)
    class(binomial_distribution), intent(in) :: this
    real(real64) :: rst
    rst = (this%n + 1.0d0) * this%p
end function

! ------------------------------------------------------------------------------
pure module function bd_variance(this) result(rst)
    class(binomial_distribution), intent(in) :: this
    real(real64) :: rst
    rst = this%n * this%p * (1.0d0 - this%p)
end function

! ------------------------------------------------------------------------------
end submodule
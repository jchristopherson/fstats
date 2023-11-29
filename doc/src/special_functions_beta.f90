submodule (fstats) special_functions_beta
contains
! ------------------------------------------------------------------------------
pure elemental module function beta_real64(a, b) result(rst)
    ! Local Variables
    real(real64), intent(in) :: a, b
    real(real64) :: rst

    ! Process
    ! REF: https://en.wikipedia.org/wiki/Beta_function
    rst = exp(log_gamma(a) + log_gamma(b) - log_gamma(a + b))
end function

! --------------------
pure elemental module function beta_real32(a, b) result(rst)
    ! Local Variables
    real(real32), intent(in) :: a, b
    real(real32) :: rst

    ! Process
    ! REF: https://en.wikipedia.org/wiki/Beta_function
    rst = exp(log_gamma(a) + log_gamma(b) - log_gamma(a + b))
end function

! ------------------------------------------------------------------------------
! source: https://people.math.sc.edu/Burkardt/f_src/special_functions/special_functions.f90
pure elemental module function regularized_beta_real64(a, b, x) result(rst)
    ! Arguments
    real(real64), intent(in) :: a, b, x
    real(real64) :: rst

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
            fk(2*k) = k * (a - k) * (one - x) / (b + two * k - one) / (b + two * k)
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

! --------------------
pure elemental module function regularized_beta_real32(a, b, x) result(rst)
    ! Arguments
    real(real32), intent(in) :: a, b, x
    real(real32) :: rst

    ! Local Variables
    real(real32), parameter :: zero = 0.0
    real(real32), parameter :: one = 1.0
    real(real32), parameter :: two = 2.0
    real(real32) :: bt, dk(51), fk(51), s0, t1, t2, ta, tb
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
            fk(2*k) = k * (a - k) * (one - x) / (b + two * k - one) / (b + two * k)
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
pure elemental module function incomplete_beta_real64(a, b, x) result(rst)
    ! Arguments
    real(real64), intent(in) :: a, b, x
    real(real64) :: rst

    ! Process
    rst = beta(a, b) * regularized_beta(a, b, x)
end function

! --------------------
pure elemental module function incomplete_beta_real32(a, b, x) result(rst)
    ! Arguments
    real(real32), intent(in) :: a, b, x
    real(real32) :: rst

    ! Process
    rst = beta(a, b) * regularized_beta(a, b, x)
end function

! ------------------------------------------------------------------------------
end submodule
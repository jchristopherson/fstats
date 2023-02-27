submodule (fstats) statistics_implementation
    use linalg, only : sort
    use ieee_arithmetic
contains
! ------------------------------------------------------------------------------
pure module function mean_real64(x) result(rst)
    ! Arguments
    real(real64), intent(in) :: x(:)
    real(real64) :: rst

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    integer(int32) :: i, n

    ! Process
    n = size(x)
    if (n == 0) then
        rst = zero
    else
        rst = x(1)
        do i = 2, n
            rst = rst + (x(i) - rst) / i
        end do
    end if
end function

! --------------------
pure module function mean_real32(x) result(rst)
    ! Arguments
    real(real32), intent(in) :: x(:)
    real(real32) :: rst

    ! Parameters
    real(real32), parameter :: zero = 0.0

    ! Local Variables
    integer(int32) :: i, n

    ! Process
    n = size(x)
    if (n == 0) then
        rst = zero
    else
        rst = x(1)
        do i = 2, n
            rst = rst + (x(i) - rst) / i
        end do
    end if
end function

! ------------------------------------------------------------------------------
pure module function variance_real64(x) result(rst)
    ! Arguments
    real(real64), intent(in) :: x(:)
    real(real64) :: rst

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: oldMean, newMean

    ! Process
    n = size(x)
    if (n <= 1) then
        rst = zero
    else
        oldMean = x(1)
        rst = zero
        do i = 2, n
            newMean = oldMean + (x(i) - oldMean) / i
            rst = rst + (x(i) - oldMean) * (x(i) - newMean)
            oldMean = newMean
        end do
        rst = rst / (n - one)
    end if
end function

! --------------------
pure module function variance_real32(x) result(rst)
    ! Arguments
    real(real32), intent(in) :: x(:)
    real(real32) :: rst

    ! Parameters
    real(real32), parameter :: zero = 0.0
    real(real32), parameter :: one = 1.0

    ! Local Variables
    integer(int32) :: i, n
    real(real32) :: oldMean, newMean

    ! Process
    n = size(x)
    if (n <= 1) then
        rst = zero
    else
        oldMean = x(1)
        rst = zero
        do i = 2, n
            newMean = oldMean + (x(i) - oldMean) / i
            rst = rst + (x(i) - oldMean) * (x(i) - newMean)
            oldMean = newMean
        end do
        rst = rst / (n - one)
    end if
end function

! ------------------------------------------------------------------------------
pure module function standard_deviation_real64(x) result(rst)
    ! Arguments
    real(real64), intent(in) :: x(:)
    real(real64) :: rst

    ! Process
    rst = sqrt(variance(x))
end function

! --------------------
pure module function standard_deviation_real32(x) result(rst)
    ! Arguments
    real(real32), intent(in) :: x(:)
    real(real32) :: rst

    ! Process
    rst = sqrt(variance(x))
end function

! ------------------------------------------------------------------------------
module function median_real64(x) result(rst)
    ! Arguments
    real(real64), intent(inout) :: x(:)
    real(real64) :: rst

    ! Parameters
    real(real64), parameter :: half = 0.5d0

    ! Local Variables
    integer(int32) :: n, nmid, nmidp1, flag, iflag

    ! Initialization
    n = size(x)
    nmid = n / 2
    nmidp1 = nmid + 1
    iflag = n - 2 * nmid
    
    ! Sort the array in ascending order
    call sort(x, .true.)

    ! Find the median
    if (iflag == 0) then
        rst = half * (x(nmid) + x(nmidp1))
    else
        rst = x(nmidp1)
    end if
end function

! --------------------
module function median_real32(x) result(rst)
    ! Arguments
    real(real32), intent(inout) :: x(:)
    real(real32) :: rst

    ! Parameters
    real(real32), parameter :: half = 0.5

    ! Local Variables
    integer(int32) :: n, nmid, nmidp1, flag, iflag

    ! Initialization
    n = size(x)
    nmid = n / 2
    nmidp1 = nmid + 1
    iflag = n - 2 * nmid
    
    ! Sort the array in ascending order
    call r32_sort(x)

    ! Find the median
    if (iflag == 0) then
        rst = half * (x(nmid) + x(nmidp1))
    else
        rst = x(nmidp1)
    end if
end function

! ------------------------------------------------------------------------------
module function r_squared_real64(x, xm, err) result(rst)
    ! Arguments
    real(real64), intent(in) :: x(:), xm(:)
    class(errors), intent(inout), optional, target :: err
    real(real64) :: rst

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: esum, vt
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = 256) :: errmsg
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    n = size(x)
    if (size(xm) /= n) then
        write(errmsg, 100) "Expected the modeled array to be of length ", &
            n, ", but found it to be of length ", size(xm), "."
        call errmgr%report_error("r_squared_real64", trim(errmsg), &
            FS_ARRAY_SIZE_ERROR)
    end if

    ! Process
    esum = zero
    do i = 1, n
        esum = esum + (x(i) - xm(i))**2
    end do
    vt = variance(x) * (n - one)
    rst = one - esum / vt

    ! Formatting
100 format(A, I0, A, I0, A)
end function

! --------------------
module function r_squared_real32(x, xm, err) result(rst)
    ! Arguments
    real(real32), intent(in) :: x(:), xm(:)
    class(errors), intent(inout), optional, target :: err
    real(real32) :: rst

    ! Parameters
    real(real32), parameter :: zero = 0.0
    real(real32), parameter :: one = 1.0

    ! Local Variables
    integer(int32) :: i, n
    real(real32) :: esum, vt
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = 256) :: errmsg
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    n = size(x)
    if (size(xm) /= n) then
        write(errmsg, 100) "Expected the modeled array to be of length ", &
            n, ", but found it to be of length ", size(xm), "."
        call errmgr%report_error("r_squared_real32", trim(errmsg), &
            FS_ARRAY_SIZE_ERROR)
    end if

    ! Process
    esum = zero
    do i = 1, n
        esum = esum + (x(i) - xm(i))**2
    end do
    vt = variance(x) * (n - one)
    rst = one - esum / vt

    ! Formatting
100 format(A, I0, A, I0, A)
end function

! ------------------------------------------------------------------------------
module function adjusted_r_squared_real64(p, x, xm, err) result(rst)
    ! Arguments
    integer(int32), intent(in) :: p
    real(real64), intent(in) :: x(:), xm(:)
    class(errors), intent(inout), optional, target :: err
    real(real64) :: rst

    ! Local Variables
    integer(int32) :: n
    real(real64) :: r2
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Parameters
    real(real64), parameter :: one = 1.0d0
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)

    ! Process
    r2 = r_squared(x, xm, errmgr)
    if (errmgr%has_error_occurred()) return
    rst = one - (one - r2) * (n - one) / (n - p - one)
end function

! --------------------
module function adjusted_r_squared_real32(p, x, xm, err) result(rst)
    ! Arguments
    integer(int32), intent(in) :: p
    real(real32), intent(in) :: x(:), xm(:)
    class(errors), intent(inout), optional, target :: err
    real(real32) :: rst

    ! Local Variables
    integer(int32) :: n
    real(real32) :: r2
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Parameters
    real(real32), parameter :: one = 1.0
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)

    ! Process
    r2 = r_squared(x, xm, errmgr)
    if (errmgr%has_error_occurred()) return
    rst = one - (one - r2) * (n - one) / (n - p - one)
end function

! ------------------------------------------------------------------------------
! REF: https://fortranwiki.org/fortran/show/Quartiles
!
! This is the method used by Minitab
pure module function quantile_real64(x, q) result(rst)
    ! Arguments
    real(real64), intent(in) :: x(:), q
    real(real64) :: rst

    ! Parameters
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    real(real64) :: a, b, c
    integer(int32) :: n, ib

    ! Initialization
    tol = sqrt(epsilon(tol))
    n = size(x)

    ! Process
    a = (n + one) * q
    b = mod(a, one)
    c = a - b

    ib = int(c, int32)
    if ((ib + 1) > n) then
        rst = (one - b) * x(ib) + b * x(n)
    else
        rst = (one - b) * x(ib) + b * x(ib + 1)
    end if
end function

! --------------------
pure module function quantile_real32(x, q) result(rst)
    ! Arguments
    real(real32), intent(in) :: x(:), q
    real(real32) :: rst

    ! Parameters
    real(real32), parameter :: one = 1.0

    ! Local Variables
    real(real32) :: a, b, c
    integer(int32) :: n, ib

    ! Initialization
    tol = sqrt(epsilon(tol))
    n = size(x)

    ! Process
    a = (n + one) * q
    b = mod(a, one)
    c = a - b

    ib = int(c, int32)
    if ((ib + 1) > n) then
        rst = (one - b) * x(ib) + b * x(n)
    else
        rst = (one - b) * x(ib) + b * x(ib + 1)
    end if
end function

! ------------------------------------------------------------------------------
pure module function confidence_interval_real64(dist, alpha, s, n) result(rst)
    ! Arguments
    class(distribution), intent(in) :: dist
    real(real64), intent(in) :: alpha, s
    integer(int32), intent(in) :: n
    real(real64) :: rst

    ! Local Variables
    integer(int32), parameter :: maxiter = 100
    real(real64), parameter :: tol = 1.0d-6
    integer(int32) :: i
    real(real64) :: x, f, df, h, twoh, dy

    ! Process
    !
    ! We use a simplified Newton's method to solve for the independent variable
    ! of the CDF function where it equals 1 - alpha / 2.
    h = 1.0d-6
    twoh = 2.0d0 * h
    x = 1.0d0 - alpha / 2.0d0
    rst = 0.5d0
    do i = 1, maxiter 
        ! Compute the CDF and its derivative at y
        f = dist%cdf(rst) - x
        df = (dist%cdf(rst + h) - dist%cdf(rst - h)) / twoh
        dy = f / df
        rst = rst - dy
        if (abs(dy) < tol) exit
    end do
end function

! ------------------
pure module function confidence_interval_real32(dist, alpha, s, n) result(rst)
    ! Arguments
    class(distribution), intent(in) :: dist
    real(real32), intent(in) :: alpha, s
    integer(int32), intent(in) :: n
    real(real32) :: rst

    ! Process
    real(real64) a, sd
    a = real(alpha, real64)
    sd = real(s, real64)
    rst = real(confidence_interval_real64(dist, a, sd, n), real32)
end function

! ------------------------------------------------------------------------------
pure module function confidence_interval_real64_array(dist, alpha, x) result(rst)
    ! Arguments
    class(distribution), intent(in) :: dist
    real(real64), intent(in) :: alpha, x(:)
    real(real64) :: rst

    ! Process
    rst = confidence_interval(dist, alpha, standard_deviation(x), size(x))
end function

! ------------------
pure module function confidence_interval_real32_array(dist, alpha, x) result(rst)
    ! Arguments
    class(distribution), intent(in) :: dist
    real(real32), intent(in) :: alpha, x(:)
    real(real32) :: rst

    ! Process
    rst = confidence_interval(dist, alpha, standard_deviation(x), size(x))
end function

! ------------------------------------------------------------------------------
module subroutine t_test_equal_var_real64(x1, x2, stat, p, dof)
    ! Arguments
    real(real64), intent(in) :: x1(:), x2(:)
    real(real64), intent(out) :: stat, p, dof

    ! Parameters
    real(real64), parameter :: half = 0.5d0
    real(real64), parameter :: one = 1.0d0
    real(real64), parameter :: two = 2.0d0

    ! Local Variables
    real(real64) :: v1, v2, m1, m2, sv, a, b, x
    integer(int32) :: n1, n2

    ! Compute the T-statistic
    n1 = size(x1)
    n2 = size(x2)
    m1 = mean(x1)
    m2 = mean(x2)
    v1 = variance(x1)
    v2 = variance(x2)
    dof = n1 + n2 - two
    sv = ((n1 - one) * v1 + (n2 - one) * v2) / dof
    stat = abs(m1 - m2) / sqrt(sv * (one / real(n1) + one / real(n2)))

    ! Compute the probability
    a = half * dof
    b = half
    x = dof / (dof + stat**2)
    p = regularized_beta(a, b, x)
end subroutine

! -------------------
module subroutine t_test_equal_var_real32(x1, x2, stat, p, dof)
    ! Arguments
    real(real32), intent(in) :: x1(:), x2(:)
    real(real32), intent(out) :: stat, p, dof

    ! Parameters
    real(real32), parameter :: half = 0.5
    real(real32), parameter :: one = 1.0
    real(real32), parameter :: two = 2.0

    ! Local Variables
    real(real32) :: v1, v2, m1, m2, sv, a, b, x
    integer(int32) :: n1, n2

    ! Compute the T-statistic
    n1 = size(x1)
    n2 = size(x2)
    m1 = mean(x1)
    m2 = mean(x2)
    v1 = variance(x1)
    v2 = variance(x2)
    dof = n1 + n2 - two
    sv = ((n1 - one) * v1 + (n2 - one) * v2) / dof
    stat = abs(m1 - m2) / sqrt(sv * (one / real(n1) + one / real(n2)))

    ! Compute the probability
    a = half * dof
    b = half
    x = dof / (dof + stat**2)
    p = regularized_beta(a, b, x)
end subroutine

! ------------------------------------------------------------------------------
module subroutine t_test_unequal_var_real64(x1, x2, stat, p, dof)
    ! Arguments
    real(real64), intent(in) :: x1(:), x2(:)
    real(real64), intent(out) :: stat, p, dof

    ! Parameters
    real(real64), parameter :: half = 0.5d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    real(real64) :: v1, v2, m1, m2, sv, a, b, x
    integer(int32) :: n1, n2

    ! Compute the T-statistic
    n1 = size(x1)
    n2 = size(x2)
    m1 = mean(x1)
    m2 = mean(x2)
    v1 = variance(x1)
    v2 = variance(x2)
    dof = (v1 / real(n1) + v2 / real(n2))**2 / ((v1 / n1)**2 / (n1 - one) + &
        (v2 / n2)**2 / (n2 - one))
    sv = sqrt(v1 / n1 + v2 / n2)
    stat = (m1 - m2) / sv

    ! Compute the probability
    a = half * dof
    b = half
    x = dof / (dof + stat**2)
    p = regularized_beta(a, b, x)
end subroutine

! -------------------
module subroutine t_test_unequal_var_real32(x1, x2, stat, p, dof)
    ! Arguments
    real(real32), intent(in) :: x1(:), x2(:)
    real(real32), intent(out) :: stat, p, dof

    ! Parameters
    real(real32), parameter :: half = 0.5
    real(real32), parameter :: one = 1.0

    ! Local Variables
    real(real32) :: v1, v2, m1, m2, sv, a, b, x
    integer(int32) :: n1, n2

    ! Compute the T-statistic
    n1 = size(x1)
    n2 = size(x2)
    m1 = mean(x1)
    m2 = mean(x2)
    v1 = variance(x1)
    v2 = variance(x2)
    dof = (v1 / real(n1) + v2 / real(n2))**2 / ((v1 / n1)**2 / (n1 - one) + &
        (v2 / n2)**2 / (n2 - one))
    sv = sqrt(v1 / n1 + v2 / n2)
    stat = (m1 - m2) / sv

    ! Compute the probability
    a = half * dof
    b = half
    x = dof / (dof + stat**2)
    p = regularized_beta(a, b, x)
end subroutine

! ------------------------------------------------------------------------------
module subroutine t_test_paired_real64(x1, x2, stat, p, dof, err)
    ! Arguments
    real(real64), intent(in) :: x1(:), x2(:)
    real(real64), intent(out) :: stat, p, dof
    class(errors), intent(inout), optional, target :: err

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: half = 0.5d0
    real(real64), parameter :: one = 1.0d0
    real(real64), parameter :: two = 2.0d0

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = 256) :: errmsg
    real(real64) :: v1, v2, m1, m2, sd, cov, a, b, x
    integer(int32) :: n1, n2, n
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n1 = size(x1)
    n2 = size(x2)
    n = min(n1, n2)

    ! Input Checking
    if (n1 /= n2) then
        write(errmsg, 100) &
            "Both arrays must be the same size.  Array 1 contains ", n1, &
            " elements and array 2 contains ", n2, "."
        call errmgr%report_error("t_test_paired_real64", trim(errmsg), &
            ML_ARRAY_SIZE_ERROR)
        return
    end if

    ! Compute the T-statistic
    m1 = mean(x1)
    m2 = mean(x2)
    v1 = variance(x1)
    v2 = variance(x2)
    dof = real(n1) - one
    cov = zero
    do i = 1, n
        cov = cov + (x1(i) - m1) * (x2(i) - m2)
    end do
    cov = cov / dof
    sd = sqrt((v1 + v2 - two * cov) / n)
    stat = (m1 - m2) / sd

    ! Compute the probability
    a = half * dof
    b = half
    x = dof / (dof + stat**2)
    p = regularized_beta(a, b, x)

    ! Formatting
100 format(A, I0, A, I0, A)
end subroutine

! -------------------
module subroutine t_test_paired_real32(x1, x2, stat, p, dof, err)
    ! Arguments
    real(real32), intent(in) :: x1(:), x2(:)
    real(real32), intent(out) :: stat, p, dof
    class(errors), intent(inout), optional, target :: err

    ! Parameters
    real(real32), parameter :: zero = 0.0
    real(real32), parameter :: half = 0.5
    real(real32), parameter :: one = 1.0
    real(real32), parameter :: two = 2.0

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = 256) :: errmsg
    real(real32) :: v1, v2, m1, m2, sd, cov, a, b, x
    integer(int32) :: n1, n2, n
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n1 = size(x1)
    n2 = size(x2)
    n = min(n1, n2)

    ! Input Checking
    if (n1 /= n2) then
        write(errmsg, 100) &
            "Both arrays must be the same size.  Array 1 contains ", n1, &
            " elements and array 2 contains ", n2, "."
        call errmgr%report_error("t_test_paired_real32", trim(errmsg), &
            ML_ARRAY_SIZE_ERROR)
        return
    end if

    ! Compute the T-statistic
    m1 = mean(x1)
    m2 = mean(x2)
    v1 = variance(x1)
    v2 = variance(x2)
    dof = real(n1) - one
    cov = zero
    do i = 1, n
        cov = cov + (x1(i) - m1) * (x2(i) - m2)
    end do
    cov = cov / dof
    sd = sqrt((v1 + v2 - two * cov) / n)
    stat = (m1 - m2) / sd

    ! Compute the probability
    a = half * dof
    b = half
    x = dof / (dof + stat**2)
    p = regularized_beta(a, b, x)

    ! Formatting
100 format(A, I0, A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
module subroutine f_test_real64(x1, x2, stat, p, dof1, dof2)
    ! Arguments
    real(real64), intent(in) :: x1(:), x2(:)
    real(real64), intent(out) :: stat, p, dof1, dof2

    ! Parameters
    real(real64), parameter :: half = 0.5d0
    real(real64), parameter :: one = 1.0d0
    real(real64), parameter :: two = 2.0d0

    ! Local Variables
    integer(int32) :: n1, n2
    real(real64) :: v1, v2, m1, m2, a, b, x

    ! Compute the F-statistic
    n1 = size(x1)
    n2 = size(x2)
    m1 = mean(x1)
    m2 = mean(x2)
    v1 = variance(x1)
    v2 = variance(x2)
    if (v1 > v2) then
        stat = v1 / v2
        dof1 = n1 - one
        dof2 = n2 - one
    else
        stat = v2 / v1
        dof1 = n2 - one
        dof2 = n1 - one
    end if

    ! Compute the probability
    a = half * dof2
    b = half * dof1
    x = dof2 / (dof2 + dof1 * stat)
    p = two * regularized_beta(a, b, x)
    if (p > one) p = two - p
end subroutine

! -------------------
module subroutine f_test_real32(x1, x2, stat, p, dof1, dof2)
    ! Arguments
    real(real32), intent(in) :: x1(:), x2(:)
    real(real32), intent(out) :: stat, p, dof1, dof2

    ! Parameters
    real(real32), parameter :: half = 0.5
    real(real32), parameter :: one = 1.0
    real(real32), parameter :: two = 2.0

    ! Local Variables
    integer(int32) :: n1, n2
    real(real32) :: v1, v2, m1, m2, a, b, x

    ! Compute the F-statistic
    n1 = size(x1)
    n2 = size(x2)
    m1 = mean(x1)
    m2 = mean(x2)
    v1 = variance(x1)
    v2 = variance(x2)
    if (v1 > v2) then
        stat = v1 / v2
        dof1 = n1 - one
        dof2 = n2 - one
    else
        stat = v2 / v1
        dof1 = n2 - one
        dof2 = n1 - one
    end if

    ! Compute the probability
    a = half * dof2
    b = half * dof1
    x = dof2 / (dof2 + dof1 * stat)
    p = two * regularized_beta(a, b, x)
    if (p > one) p = two - p
end subroutine

! ------------------------------------------------------------------------------
! REF: https://www.spcforexcel.com/knowledge/root-cause-analysis/single-factor-anova
module function anova_1_factor(x) result(rst)
    ! Arguments
    real(real64), intent(in) :: x(:,:)
    type(single_factor_anova_table) :: rst

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    integer(int32) :: j, a, n, nt
    real(real64) :: sum_all, tssq, essq, bssq

    ! Initialization
    a = size(x, 2)
    nt = size(x, 1)
    n = nt * a
    rst%within_factor%f_statistic = ieee_value(sum_all, IEEE_QUIET_NAN)
    rst%within_factor%probability = ieee_value(sum_all, IEEE_QUIET_NAN)

    ! Determine the degrees of freedom
    rst%main_factor%dof = a - 1
    rst%within_factor%dof = n - a
    rst%total_dof = n - 1

    ! Quick Return
    if (a == 1 .or. nt == 1) then
        rst%main_factor%sum_of_squares = zero
        rst%main_factor%variance = zero
        rst%main_factor%f_statistic = zero
        rst%main_factor%probability = zero
        rst%within_factor%sum_of_squares = zero
        rst%within_factor%variance = zero
        rst%total_variance = variance(pack(x, .true.))
        rst%total_sum_of_squares = rst%total_variance * rst%total_dof
        rst%overall_mean = mean(pack(x, .true.))
        return
    end if

    ! Compute the sum of squares for all factors
    sum_all = sum(x)
    tssq = sum(x**2) - (sum_all**2 / n)
    
    bssq = zero
    do j = 1, a
        bssq = bssq + sum(x(:,j))**2
    end do
    bssq = (bssq / nt) - (sum_all**2 / n)
    essq = tssq - bssq

    rst%main_factor%sum_of_squares = bssq
    rst%within_factor%sum_of_squares = essq
    rst%total_sum_of_squares = tssq

    ! Compute the variance terms
    rst%main_factor%variance = bssq / rst%main_factor%dof
    rst%within_factor%variance = essq / rst%within_factor%dof
    rst%total_variance = tssq / rst%total_dof

    ! Compute the overall mean
    rst%overall_mean = mean(pack(x, .true.))

    ! Compute the F-statistic and probability term
    call anova_probability( &
        rst%main_factor%variance, &
        rst%within_factor%variance, &
        rst%main_factor%dof, &
        rst%within_factor%dof, &
        rst%main_factor%f_statistic, &
        rst%main_factor%probability &
    )
end function

! ------------------------------------------------------------------------------
! REF: https://www.spcforexcel.com/knowledge/measurement-systems-analysis/anova-gage-rr-part-1
! REF: https://www.itl.nist.gov/div898/handbook/prc/section4/prc427.htm
! Data set is expected as a 3D array with each of the K pages containing the R 
!   treatments of N tests such that the array size is N-by-R-by-K
module function anova_2_factor(x) result(rst)
    ! Arguments
    real(real64), intent(in) :: x(:,:,:)
    type(two_factor_anova_table) :: rst

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, j, jj, k, r, n
    real(real64) :: factorMean
    real(real64), allocatable :: xpack(:)

    ! Initialization
    n = size(x, 3)
    k = size(x, 2)
    r = size(x, 1)
    rst%within_factor%f_statistic = ieee_value(sum_all, IEEE_QUIET_NAN)
    rst%within_factor%probability = ieee_value(sum_all, IEEE_QUIET_NAN)

    ! Quick Return
    if (k == 1) then
        ! This is a one-factor anova
    end if

    ! Determine the number of DOF
    rst%main_factor_1%dof = k - one
    rst%main_factor_2%dof = n - 1
    rst%interaction%dof = (k - 1) * (n - 1)
    rst%within_factor%dof = n * k * (r - 1)
    rst%total_dof = n * k * r - 1

    ! Compute the overall mean, sum of squares, and variance
    xpack = pack(x, .true.)
    rst%overall_mean = mean(xpack)
    rst%total_sum_of_squares = sum((xpack - rst%overall_mean)**2)
    rst%total_variance = rst%total_sum_of_squares / rst%total_dof

    ! Compute factor 1 results
    rst%main_factor_1%sum_of_squares = zero
    do i = 1, k
        factorMean = mean(pack(x(:,i,:), .true.))
        rst%main_factor_1%sum_of_squares = rst%main_factor_1%sum_of_squares + &
            (factorMean - rst%overall_mean)**2
    end do
    rst%main_factor_1%sum_of_squares = n * r * rst%main_factor_1%sum_of_squares
    rst%main_factor_1%variance = rst%main_factor_1%sum_of_squares / &
        rst%main_factor_1%dof

    ! Compute factor 2 results
    rst%main_factor_2%sum_of_squares = zero
    do i = 1, n
        factorMean = mean(pack(x(:,:,i), .true.))
        rst%main_factor_2%sum_of_squares = rst%main_factor_2%sum_of_squares + &
            (factorMean - rst%overall_mean)**2
    end do
    rst%main_factor_2%sum_of_squares = k * r * rst%main_factor_2%sum_of_squares
    rst%main_factor_2%variance = rst%main_factor_2%sum_of_squares / &
        rst%main_factor_2%dof

    ! Compute the within (error) term
    rst%within_factor%sum_of_squares = zero
    do j = 1, k
        do i = 1, n
            factorMean = mean(x(:,j,i))
            do jj = 1, r
                rst%within_factor%sum_of_squares = &
                    rst%within_factor%sum_of_squares + &
                    (x(jj,j,i) - factorMean)**2
            end do
        end do
    end do
    rst%within_factor%variance = rst%within_factor%sum_of_squares /&
         rst%within_factor%dof

    ! Compute the interaction term
    rst%interaction%sum_of_squares = rst%total_sum_of_squares - ( &
        rst%main_factor_1%sum_of_squares + &
        rst%main_factor_2%sum_of_squares + &
        rst%within_factor%sum_of_squares &
    )
    rst%interaction%variance = rst%interaction%sum_of_squares / &
        rst%interaction%dof

    ! Compute the F-statistics
    call anova_probability( &
        rst%main_factor_1%variance, &
        rst%within_factor%variance, &
        rst%main_factor_1%dof, &
        rst%within_factor%dof, &
        rst%main_factor_1%f_statistic, &
        rst%main_factor_1%probability &
    )
    call anova_probability( &
        rst%main_factor_2%variance, &
        rst%within_factor%variance, &
        rst%main_factor_2%dof, &
        rst%within_factor%dof, &
        rst%main_factor_2%f_statistic, &
        rst%main_factor_2%probability &
    )
    call anova_probability( &
        rst%interaction%variance, &
        rst%within_factor%variance, &
        rst%interaction%dof, &
        rst%within_factor%dof, &
        rst%interaction%f_statistic, &
        rst%interaction%probability &
    )
end function

! ------------------------------------------------------------------------------
! REF: https://www.spcforexcel.com/knowledge/root-cause-analysis/understanding-regression-statistics-part-1
module function anova_model_fit(nmodelparams, ymeas, ymod, err) result(rst)
    ! Arguments
    integer(int32), intent(in) :: nmodelparams
    real(real64), intent(in) :: ymeas(:), ymod(:)
    class(errors), intent(inout), optional, target :: err
    type(single_factor_anova_table) :: rst

    ! Local Variables
    integer(int32) :: n, flag
    real(real64), allocatable :: ypack(:)
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = 256) :: errmsg
    
    ! Initialization
    n = size(ymeas)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    rst%within_factor%f_statistic = ieee_value(sum_all, IEEE_QUIET_NAN)
    rst%within_factor%probability = ieee_value(sum_all, IEEE_QUIET_NAN)

    ! Input Checking
    if (size(ymod) /= n) then
        write(errmsg, 100) &
            "Both input arrays must be the same size.  Expected ", n, &
            ", but found ", size(ymod), "."
        call errmgr%report_error("anova_model_fit", trim(errmsg), &
            ML_ARRAY_SIZE_ERROR)
        return
    end if

    ! Memory Allocation
    allocate(ypack(2 * n), stat = flag)
    if (flag /= 0) then
        write(errmsg, 101) "Memory allocation error code ", flag, "."
        call errmgr%report_error("anova_model_fit", &
            trim(errmsg), ML_OUT_OF_MEMORY_ERROR)
        return
    end if

    ! Determine the number of DOF
    rst%main_factor%dof = nmodelparams - 1
    rst%within_factor%dof = n - rst%main_factor%dof - 1
    rst%total_dof = n - 1

    ! Process
    ypack(1:n) = ymeas
    ypack(n+1:2*n) = ymod
    rst%overall_mean = mean(ypack)
    rst%total_sum_of_squares = sum((ymeas - rst%overall_mean)**2)
    rst%main_factor%sum_of_squares = sum((ymod - rst%overall_mean)**2)
    rst%within_factor%sum_of_squares = sum((ymeas - ymod)**2)

    rst%total_variance = rst%total_sum_of_squares / rst%total_dof
    rst%main_factor%variance = rst%main_factor%sum_of_squares / &
        rst%main_factor%dof
    rst%within_factor%variance = rst%within_factor%sum_of_squares / &
        rst%within_factor%dof

    ! Compute the F-statistic and probability term
        call anova_probability( &
        rst%main_factor%variance, &
        rst%within_factor%variance, &
        rst%main_factor%dof, &
        rst%within_factor%dof, &
        rst%main_factor%f_statistic, &
        rst%main_factor%probability &
    )    

    ! Formatting
100 format(A, I0, A, I0, A)
101 format(A, I0, A)
end function

! ******************************************************************************
! PRIVATE ROUTINES
! ------------------------------------------------------------------------------
! This is a hack as the linalg library doesn't support 32-bit real32ing point
! routines
subroutine r32_sort(x)
    ! Arguments
    real(real32), intent(inout) :: x(:)

    ! Local Variables
    real(real64), allocatable :: xc(:)

    ! Process
    xc = real(x, real64)
    call sort(xc, .true.)
    x = real(xc, real32)
end subroutine

! ------------------------------------------------------------------------------
subroutine anova_probability(v1, v2, dof1, dof2, f, p)
    ! Arguments
    real(real64), intent(in) :: v1, v2, dof1, dof2
    real(real64), intent(out) :: f, p

    ! Local Variables
    real(real64) :: d1, d2, a, b, x
    
    ! Process
    f = v1 / v2
    d1 = dof1
    d2 = dof2

    a = 0.5d0 * d2
    b = 0.5d0 * d1
    x = d2 / (d2 + d1 * f)

    p = regularized_beta(a, b, x)
    if (p > 1.0d0) then
        p = 2.0d0 - p
    end if
end subroutine

! ------------------------------------------------------------------------------
end submodule
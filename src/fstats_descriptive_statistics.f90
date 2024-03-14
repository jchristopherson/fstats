module fstats_descriptive_statistics
    use linalg, only : sort
    use ieee_arithmetic
    use fstats_errors
    use fstats_special_functions
    use fstats_distributions
    use ferror
    implicit none
    private
    public :: mean
    public :: variance
    public :: standard_deviation
    public :: median
    public :: r_squared
    public :: adjusted_r_squared
    public :: quantile
    public :: confidence_interval
    public :: t_test_equal_variance
    public :: t_test_unequal_variance
    public :: t_test_paired
    public :: f_test
    public :: trimmed_mean

    interface confidence_interval
        !! Computes the confidence interval for the specified distribution.
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Confidence_interval)
        module procedure :: confidence_interval_scalar
        module procedure :: confidence_interval_array
    end interface
contains
! ------------------------------------------------------------------------------
pure function mean(x) result(rst)
    !! Computes the mean of the values in an array.
    real(real64), intent(in) :: x(:)
        !! The array of values to analyze.
    real(real64) :: rst
        !! The result.

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

! ------------------------------------------------------------------------------
pure function variance(x) result(rst)
    !! Computes the sample variance of the values in an array.
    !!
    !! The variance computed is the sample variance such that the variance.
    !! $$ s^2 = \frac{\Sigma \left( x_{i} - \bar{x} \right)^2}{n - 1} $$
    real(real64), intent(in) :: x(:)
        !! The array of values to analyze.
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

! ------------------------------------------------------------------------------
pure function standard_deviation(x) result(rst)
    !! Computes the sample standard deviation of the values in an array.
    !! 
    !! The value computed is the sample standard deviation.
    !! $$ s = \sqrt{ \frac{\Sigma \left( x_{i} - \bar{x} \right)^2}{n - 1} } $$
    real(real64), intent(in) :: x(:)
        !! The array of values to analyze.
    real(real64) :: rst
        !! The result.

    ! Process
    rst = sqrt(variance(x))
end function

! ------------------------------------------------------------------------------
function median(x) result(rst)
    !! Computes the median of the values in an array.
    real(real64), intent(inout) :: x(:)
        !! The array of values to analyze.  On output, this array is sorted into
        !! ascending order.
    real(real64) :: rst
        !! The result.

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


! ------------------------------------------------------------------------------
function r_squared(x, xm, err) result(rst)
    !! Computes the R-squared value for a data set.
    !!
    !! The R-squared value is computed by determining the sum of the squares
    !! of the residuals: 
    !! $$ SS_{res} = \Sigma \left( y_i - f_i \right)^2 $$
    !! The total sum of the squares: 
    !! $$ SS_{tot} = \Sigma \left( y_i - \bar{y} \right)^2 $$. 
    !! The R-squared value is then: 
    !! $$ R^2 = 1 - \frac{SS_{res}}{SS_{tot}} $$.
    !!
    !! See Also:
    !!
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Coefficient_of_determination)
    real(real64), intent(in) :: x(:)
        !! An N-element array containing the dependent variables from 
        !! the data set.
    real(real64), intent(in) :: xm(:)
        !! An N-element array containing the corresponding modeled 
        !! values.
    class(errors), intent(inout), optional, target :: err
        !! A mechanism for communicating errors and warnings
        !! to the caller.  Possible warning and error codes are as 
        !! follows.
        !! - FS_NO_ERROR: No errors encountered.
        !! - FS_ARRAY_SIZE_ERROR: Occurs if x and xm are not the 
        !!   same size.
    real(real64) :: rst
        !! The result.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: esum, vt
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    n = size(x)
    if (size(xm) /= n) then
        call report_array_size_error(errmgr, "r_squared_real64", "XM", n, &
            size(xm))
        return
    end if

    ! Process
    esum = zero
    do i = 1, n
        esum = esum + (x(i) - xm(i))**2
    end do
    vt = variance(x) * (n - one)
    rst = one - esum / vt
end function

! ------------------------------------------------------------------------------
function adjusted_r_squared(p, x, xm, err) result(rst)
    !! Computes the adjusted R-squared value for a data set.
    !!
    !! The adjusted R-squared provides a mechanism for tempering the effects
    !! of extra explanatory variables on the traditional R-squared 
    !! calculation.  It is computed by noting the sample size \( n \) and 
    !! the number of variables \( p \).
    !! $$ \bar{R}^2 = 1 - \left( 1 - R^2 \right) \frac{n - 1}{n - p} $$.
    !!
    !! See Also:
    !!
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Coefficient_of_determination#Adjusted_R2)
    integer(int32), intent(in) :: p
        !! The number of variables.
    real(real64), intent(in) :: x(:)
        !! An N-element array containing the dependent variables from 
        !! the data set.
    real(real64), intent(in) :: xm(:)
        !! An N-element array containing the corresponding modeled 
        !! values.
    class(errors), intent(inout), optional, target :: err
        !! A mechanism for communicating errors and warnings
        !! to the caller.  Possible warning and error codes are as 
        !! follows.
        !! - FS_NO_ERROR: No errors encountered.
        !! - FS_ARRAY_SIZE_ERROR: Occurs if x and xm are not the 
        !!   same size.
    real(real64) :: rst
        !! The result.

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

! ------------------------------------------------------------------------------
! REF: https://fortranwiki.org/fortran/show/Quartiles
!
! This is the method used by Minitab
pure function quantile(x, q) result(rst)
    !! Computes the specified quantile of a data set using the SAS 
    !! Method 4.
    !!
    !! See Also
    !!
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Quantile)
    real(real64), intent(in) :: x(:)
        !! An N-element array containing the data.
    real(real64), intent(in) :: q
        !! The quantile to compute (e.g. 0.25 computes the 25% quantile).
    real(real64) :: rst
        !! The result.

    ! Parameters
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    real(real64) :: a, b, c, tol
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
pure function confidence_interval_scalar(dist, alpha, s, n) result(rst)
    !! Computes the confidence interval for the specified distribution.
    class(distribution), intent(in) :: dist
        !! The distribution object defining the probability distribution
        !! to establish the confidence level.
    real(real64), intent(in) :: alpha
        !! The probability value of interest.  For instance, use a value of 0.05
        !! for a confidence level of 95%.
    real(real64), intent(in) :: s
        !! The sample standard deviation.
    integer(int32), intent(in) :: n
        !! The number of samples in the data set.
    real(real64) :: rst
        !! The result.

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

    ! Determine the actual interval
    rst = rst * s / sqrt(real(n, real64))
end function

! ------------------------------------------------------------------------------
pure function confidence_interval_array(dist, alpha, x) result(rst)
    !! Computes the confidence interval for the specified distribution.
    class(distribution), intent(in) :: dist
        !! The distribution object defining the probability distribution
        !! to establish the confidence level.
    real(real64), intent(in) :: alpha
        !! The probability value of interest.  For instance, use a value of 0.05
        !! for a confidence level of 95%.
    real(real64), intent(in) :: x(:)
        !! An N-element array containing the data to analyze.
    real(real64) :: rst
        !! The result.

    ! Process
    rst = confidence_interval(dist, alpha, standard_deviation(x), size(x))
end function

! ------------------------------------------------------------------------------
subroutine t_test_equal_variance(x1, x2, stat, p, dof)
    !! Computes the 2-tailed Student's T-Test for two data sets of 
    !! assumed equivalent variances.
    !!
    !! See Also
    !!
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Student%27s_t-test)
    real(real64), intent(in) :: x1(:)
        !! An N-element array containing the first data set.
    real(real64), intent(in) :: x2(:)
        !! An M-element array containing the second data set.
    real(real64), intent(out) :: stat
        !! The Student-'s T-Test statistic.
    real(real64), intent(out) :: p
        !! The probability value that the two samples are likely to
        !! have come from two underlying populations that 
        !! have the same mean.
    real(real64), intent(out) :: dof
        !! The degrees of freedom.

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

! ------------------------------------------------------------------------------
subroutine t_test_unequal_variance(x1, x2, stat, p, dof)
    !! Computes the 2-tailed Student's T-Test for two data sets of 
    !! assumed non-equivalent variances.
    !!
    !! See Also
    !!
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Student%27s_t-test)
    real(real64), intent(in) :: x1(:)
        !! An N-element array containing the first data set.
    real(real64), intent(in) :: x2(:)
        !! An M-element array containing the second data set.
    real(real64), intent(out) :: stat
        !! The Student-'s T-Test statistic.
    real(real64), intent(out) :: p
        !! The probability value that the two samples are likely to
        !! have come from two underlying populations that 
        !! have the same mean.
    real(real64), intent(out) :: dof
        !! The degrees of freedom.

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

! ------------------------------------------------------------------------------
subroutine t_test_paired(x1, x2, stat, p, dof, err)
    !! Computes the 2-tailed Student's T-Test for two paired data sets.
    !!
    !! See Also
    !!
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Student%27s_t-test)
    real(real64), intent(in) :: x1(:)
        !! An N-element array containing the first data set.
    real(real64), intent(in) :: x2(:)
        !! An N-element array containing the second data set.
    real(real64), intent(out) :: stat
        !! The Student-'s T-Test statistic.
    real(real64), intent(out) :: p
        !! The probability value that the two samples are likely to
        !! have come from  two underlying populations that 
        !! have the same mean.
    real(real64), intent(out) :: dof
        !! The degrees of freedom.
    class(errors), intent(inout), optional, target :: err
        !! A mechanism for communicating errors and warnings to the 
        !! caller.  Possible warning and error codes are as follows.
        !! - FS_NO_ERROR: No errors encountered.
        !! - FS_ARRAY_SIZE_ERROR: Occurs if x1 and x2 are not the same 
        !!   length.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: half = 0.5d0
    real(real64), parameter :: one = 1.0d0
    real(real64), parameter :: two = 2.0d0

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    real(real64) :: v1, v2, m1, m2, sd, cov, a, b, x
    integer(int32) :: i, n1, n2, n
    
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
        call report_arrays_not_same_size_error(errmgr, "t_test_paired_real64", &
            "X1", "X2", n1, n2)
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
end subroutine

! ------------------------------------------------------------------------------
subroutine f_test(x1, x2, stat, p, dof1, dof2)
    !! Computes the F-test and returns the probability (two-tailed) that
    !! the variances of two data sets are not significantly different.
    !!
    !! See Also
    !!
    !! - [Wikipedia](https://en.wikipedia.org/wiki/F-test)
    real(real64), intent(in) :: x1(:)
        !! An N-element array containing the first data set.
    real(real64), intent(in) :: x2(:)
        !! An M-element array containing the second data set.
    real(real64), intent(out) :: stat
        !! The F-statistic.
    real(real64), intent(out) :: p
        !! The probability value that the two samples are likely to
        !! have come from the two underlying populations that 
        !! have the same variance.
    real(real64), intent(out) :: dof1
        !! A measure of the degrees of freedom.
    real(real64), intent(out) :: dof2
        !! A measure of the degrees of freedom.

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

! ------------------------------------------------------------------------------
function trimmed_mean(x, p) result(rst)
    !! Computes the trimmed mean of a data set.
    real(real64), intent(inout), dimension(:) :: x
        !! An N-element array containing the data.  On output, the
        !! array is sorted into ascending order.
    real(real64), intent(in), optional :: p
        !! An optional parameter specifying the percentage of values
        !! from either end of the distribution to remove.  The default
        !! is 0.05 such that the bottom 5% and top 5% are removed.
    real(real64) :: rst
        !! The trimmed mean.

    ! Local Variables
    integer(int32) :: i1, i2, n
    real(real64) :: pv

    ! Initialization
    if (present(p)) then
        pv = abs(p)
    else
        pv = 0.05d0
    end if

    ! Sort the array into ascending order
    call sort(x, .true.)

    ! Find the limiting indices
    n = size(x)
    i1 = max(floor(n * pv, int32), 1)
    i2 = min(n, n - i1 + 1)
    rst = mean(x(i1:i2))
end function

! ------------------------------------------------------------------------------
end module
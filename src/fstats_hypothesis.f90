module fstats_hypothesis
    use iso_fortran_env
    use ieee_arithmetic
    use fstats_errors
    use fstats_special_functions
    use fstats_distributions
    use fstats_descriptive_statistics
    use fstats_types
    private
    public :: confidence_interval
    public :: t_test_equal_variance
    public :: t_test_unequal_variance
    public :: t_test_paired
    public :: f_test
    public :: bartletts_test
    public :: levenes_test
    public :: sample_size

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
    real(real64) :: x, dn

    ! Process
    dn = real(n, real64)
    x = 1.0d0 - 0.5d0 * alpha
    rst = dist%standardized_variable(x)
    rst = rst * s / sqrt(dn)
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
    real(real64), parameter :: one = 1.0d0
    real(real64), parameter :: two = 2.0d0

    ! Local Variables
    integer(int32) :: n1, n2
    real(real64) :: v1, v2, m1, m2
    type(f_distribution) :: dist

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

    dist%d1 = dof1
    dist%d2 = dof2
    p = two * (one - dist%cdf(stat))! 2x because this is a two-tailed estimate
    if (p > one) p = two - p
end subroutine

! ------------------------------------------------------------------------------
subroutine bartletts_test(x, stat, p)
    !! Computes Bartlett's test statistic and associated probability.
    !!
    !! The statistic is calculated as follows.
    !!
    !! $$ \chi^{2} = \frac{(N - k) \ln(S_{p}^{2}) \sum_{i = 1}^{k} 
    !! \left(n_{i} - 1 \right) \ln(S_{i}^{2})}{1 + 
    !! \frac{1}{3 \left( k - 1 \right)} \left( \sum_{i = 1}^{k} 
    !! \left( \frac{1}{n_{i} - 1} \right) - \frac{1}{N - k} \right)} $$
    !!
    !! Where \( N = \sum_{i = 1}^{k} n_{i} \) and \( S_{p}^{2} \) is the pooled
    !! variance.
    !!
    !! The probability is calculated as the right-tail probability of the
    !! chi-squared distribution.
    !!
    !! Bartlett's test is most relevant for distributions showing strong 
    !! normality.  For distributions lacking strong normality, consider 
    !! Levene's test instead.
    !!
    !! See Also
    !!
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Bartlett%27s_test)
    type(array_container), intent(in), dimension(:) :: x
        !! The arrays of data to analyze.
    real(real64), intent(out) :: stat
        !! The Bartlett's test statistic.
    real(real64), intent(out) :: p
        !! The probability value that the variances of each data set are
        !! equivalent.  A low p-value, less than some significance level,
        !! indicates a non-equivalance of variances.

    ! Local Variables
    integer(int32) :: i, n, k, ni
    real(real64) :: si, sp, numer, denom
    type(chi_squared_distribution) :: dist

    ! Initialization
    k = size(x)
    n = 0
    do i = 1, k
        n = n + size(x(i)%x)
    end do

    ! Compute the statistic
    n = 0
    sp = 0.0d0
    numer = 0.0d0
    denom = 0.0d0
    do i = 1, k
        ni = size(x(i)%x)
        n = n + ni
        si = variance(x(i)%x)
        sp = sp + (ni - 1.0d0) * si
        numer = numer + (ni - 1.0d0) * log(variance(x(i)%x))
        denom = denom + 1.0d0 / (ni - 1.0d0)
    end do
    sp = sp / real(n - k, real64)
    stat = ((n - k) * log(sp) - numer) / &
        (1.0d0 + (1.0d0 / (3.0d0 * k - 3.0d0)) * &
        (denom - 1.0d0 / real(n - k, real64)))

    ! Compute the p-value
    dist%dof = k - 1
    p = 1.0d0 - dist%cdf(stat)
end subroutine

! ------------------------------------------------------------------------------
subroutine levenes_test(x, stat, p, err)
    !! Computes Levene's test statistic and associated probability.
    !!
    !! The statistic is calculated as follows.
    !! $$ W = \frac{N - k}{k - 1} \frac{ \sum_{i = 1}^{k} N_{i} \left( Z_{i.} - 
    !! Z{..} \right)^{2}}{ \sum_{i = 1}^{k} \sum_{j = 1}^{n_{i}} \left( Z_{ij} -
    !! Z_{i.} \right)^{2} } $$
    !!
    !! Where:
    !! $$ Z_{ij} = |X_{ij} - \overline{X_{i.}}| $$
    !! $$ Z_{i.} = \frac{1}{n_{i}} \sum_{j = 1}^{n_{i}} Z_{ij} $$ 
    !! $$ Z_{..} = \frac{1}{N} \sum_{i = 1}^{k} \sum_{j = 1}^{n_{i}} Z_{ij} $$
    !!
    !! As the test statistic is approximately F-distributed, the F-distribution
    !! is used to calculate the probability term.
    !!
    !! See Also
    !!
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Levene%27s_test)
    type(array_container), intent(in), dimension(:) :: x
        !! The arrays of data to analyze.
    real(real64), intent(out) :: stat
        !! The Bartlett's test statistic.
    real(real64), intent(out) :: p
        !! The probability value that the variances of each data set are
        !! equivalent.  A low p-value, less than some significance level,
        !! indicates a non-equivalance of variances.
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    integer(int32) :: i, j, k, n, ni, flag
    real(real64) :: numer, denom, inner, yi, z, zij
    real(real64), allocatable, dimension(:) :: y, zt, zi
    type(f_distribution) :: dist
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    k = size(x)

    ! Local Memory Allocations
    allocate(y(k), zi(k), stat = flag)
    if (flag /= 0) then
        call report_memory_error(errmgr, "levenes_test", flag)
        return
    end if

    ! Compute the total mean
    z = 0.0d0
    n = 0
    do i = 1, k
        ni = size(x(i)%x)
        n = n + ni
        y(i) = mean(x(i)%x)
        zt = abs(x(i)%x - y(i))
        zi(i) = mean(zt)
        z = z + zi(i) * ni
    end do
    z = z / n

    ! Process
    numer = 0.0d0
    denom = 0.0d0
    do i = 1, k
        ni = size(x(i)%x)
        yi = y(i)
        numer = numer + ni * (zi(i) - z)**2
        
        inner = 0.0d0
        do j = 1, ni
            zij = abs(x(i)%x(j) - yi)
            inner = inner + (zij - zi(i))**2
        end do
        denom = denom + inner
    end do
    stat = real((N - k) / (k - 1), real64) * (numer / denom)
    dist%d1 = k - 1.0d0
    dist%d2 = real(n - k, real64)
    p = 1.0d0 - dist%cdf(stat)
end subroutine

! ------------------------------------------------------------------------------
pure function sample_size(dist, var, delta, bet, alpha) result(rst)
    !! Estimates the sample size required to achieve an experiment with the
    !! desired power and significance levels to ascertain the desired 
    !! difference in parameter.
    !!
    !! See Also
    !!
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Power_of_a_test)
    class(distribution), intent(in) :: dist
        !! The distribution to utilize as a measure.
    real(real64), intent(in) :: var
        !! An estimate of the population variance.
    real(real64), intent(in) :: delta
        !! The parameter difference that is desired.
    real(real64), intent(in), optional :: bet
        !! The desired power level.  The default for this value is 0.2, for a 
        !! power of 80%.
    real(real64), intent(in), optional :: alpha
        !! The desired significance level.  The default for this value is 0.05
        !! for a confidence level of 95%.
    real(real64) :: rst
        !! The minimum sample size requried to achieve the desired experimental
        !! outcome.

    ! Local Variables
    real(real64) :: a, b, za, zb

    ! Initialization
    if (present(bet)) then
        b = bet
    else
        b = 0.8d0
    end if
    if (present(alpha)) then
        a = alpha
    else
        a = 0.05d0
    end if

    za = dist%standardized_variable(1.0d0 - a / 2.0d0)
    zb = dist%standardized_variable(b)
    rst = 2.0d0 * (za + zb)**2 * var / (delta**2)
end function

! ------------------------------------------------------------------------------
end module
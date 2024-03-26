module fstats_descriptive_statistics
    use iso_fortran_env
    use linalg, only : sort
    use ferror
    use fstats_errors
    implicit none
    private
    public :: mean
    public :: variance
    public :: standard_deviation
    public :: median
    public :: quantile
    public :: trimmed_mean
    public :: covariance
    
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
    !! The variance computed is the sample variance such that 
    !! $$ s^2 = \frac{\Sigma \left( x_{i} - \bar{x} \right)^2}{n - 1} $$.
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
pure function covariance(x, y) result(rst)
    !! Computes the sample covariance of two data sets.
    !!
    !! The covariance computed is the sample covariance such that 
    !! $$ q_{jk} = \frac{\Sigma \left( x_{i} - \bar{x} \right) 
    !! \left( y_{i} - \bar{y} \right)}{n - 1} $$.
    real(real64), intent(in), dimension(:) :: x
        !! The first N-element data set.
    real(real64), intent(in), dimension(size(x)) :: y
        !! The second N-element data set.
    real(real64) :: rst
        !! The covariance.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: meanX, meanY

    ! Process
    n = size(x)
    if (n <= 1) then
        rst = zero
    else
        ! Compute the means
        meanX = x(1)
        meanY = y(1)
        do i = 2, n
            meanX = meanX + (x(i) - meanX) / i
            meanY = meanY + (y(i) - meanY) / i
        end do

        ! Compute the covariance
        rst = sum((x - meanX) * (y - meanY)) / (n - one)
    end if
end function

! ------------------------------------------------------------------------------
end module
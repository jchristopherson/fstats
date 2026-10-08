! SPDX-FileCopyrightText: 2022-2026 Jason Christopherson
! SPDX-License-Identifier: MIT
!
! MIT License
!
! Copyright (c) 2022-2026 Jason Christopherson
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module fstats_descriptive_statistics
    use iso_fortran_env
    use ieee_arithmetic
    use linalg, only : sort
    use fstats_errors
    use fstats_types
    implicit none
    private
    public :: mean
    public :: variance
    public :: standard_deviation
    public :: median
    public :: quantile
    public :: trimmed_mean
    public :: covariance
    public :: pooled_variance
    
    interface pooled_variance
        !! Computes the pooled estimate of variance.
        module procedure :: pooled_variance_1
        module procedure :: pooled_variance_2
    end interface
contains
pure function shifted_scaled_values(x, scale_value) result(rst)
    !! Shifts finite observations by their first value and scales their differences.
    !! Opposite-sign subtraction is performed after scaling to avoid overflow.
    real(real64), intent(in), dimension(:) :: x
        !! Nonempty, finite observations.
    real(real64), intent(in) :: scale_value
        !! Positive finite scale, normally MAXVAL(ABS(x)).
    real(real64), allocatable, dimension(:) :: rst
        !! Shifted/scaled observations; the first entry is zero.
    integer(int32) :: observation

    allocate(rst(size(x)))
    do observation = 1, size(x)
        if ((x(observation) >= 0.0d0 .and. x(1) >= 0.0d0) .or. &
            (x(observation) <= 0.0d0 .and. x(1) <= 0.0d0)) then
            rst(observation) = (x(observation) - x(1)) / scale_value
        else
            rst(observation) = x(observation) / scale_value - x(1) / scale_value
        end if
    end do
end function

pure function sample_spread(x) result(rst)
    !! Computes the sample standard deviation with shifted, scaled Welford updates.
    real(real64), intent(in), dimension(:) :: x
        !! Finite observations; at least two are required.
    real(real64) :: rst
        !! Sample standard deviation; invalid/undersized data give NaN,
        !! and an unrepresentable result gives positive infinity.
    real(real64), allocatable, dimension(:) :: shifted
    real(real64) :: scale_value, running_mean, delta, second_moment, factor
    integer(int32) :: observation, count

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    count = size(x)
    if (count < 2 .or. .not.all(ieee_is_finite(x))) return
    scale_value = maxval(abs(x))
    rst = 0.0d0
    if (scale_value == 0.0d0) return
    shifted = shifted_scaled_values(x, scale_value)
    running_mean = 0.0d0
    second_moment = 0.0d0
    do observation = 2, count
        delta = shifted(observation) - running_mean
        running_mean = running_mean + delta / real(observation, real64)
        second_moment = second_moment + delta * (shifted(observation) - running_mean)
    end do
    factor = sqrt(max(0.0d0, second_moment / real(count - 1, real64)))
    if (factor > 1.0d0) then
        if (scale_value > huge(1.0d0) / factor) then
            rst = ieee_value(0.0d0, ieee_positive_inf)
            return
        end if
    end if
    rst = scale_value * factor
end function

! ------------------------------------------------------------------------------
pure function mean(x) result(rst)
    !! Computes the mean using a scaled, compensated sum.
    real(real64), intent(in) :: x(:)
        !! The array of values to analyze.
    real(real64) :: rst
        !! Mean; empty arrays or nonfinite observations give NaN.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: scale_value, compensation, term, updated

    ! Process
    n = size(x)
    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (n == 0 .or. .not.all(ieee_is_finite(x))) return
    scale_value = maxval(abs(x))
    rst = zero
    if (scale_value == zero) return
    compensation = zero
    do i = 1, n
        term = x(i) / scale_value
        updated = rst + term
        if (abs(rst) >= abs(term)) then
            compensation = compensation + ((rst - updated) + term)
        else
            compensation = compensation + ((term - updated) + rst)
        end if
        rst = updated
    end do
    rst = scale_value * ((rst + compensation) / real(n, real64))
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
        !! Sample variance; fewer than two observations or nonfinite inputs
        !! give NaN. An unrepresentable variance gives positive infinity.

    rst = sample_spread(x)
    if (rst > sqrt(huge(1.0d0))) then
        rst = ieee_value(0.0d0, ieee_positive_inf)
    else
        rst = rst**2
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
        !! Sample standard deviation; invalid/undersized data give NaN.
        !! It can remain finite even when the variance would overflow.

    ! Process
    rst = sample_spread(x)
end function

! ------------------------------------------------------------------------------
pure function median(x) result(rst)
    !! Computes the median of the values in an array.
    real(real64), intent(in) :: x(:)
        !! The array of values to analyze.
    real(real64) :: rst
        !! Median; empty arrays or nonfinite observations give NaN.

    ! Parameters
    real(real64), parameter :: half = 0.5d0

    ! Local Variables
    integer(int32) :: n, nmid, nmidp1, flag, iflag
    real(real64), allocatable, dimension(:) :: xc

    ! Initialization
    n = size(x)
    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (n == 0 .or. .not.all(ieee_is_finite(x))) return
    nmid = n / 2
    nmidp1 = nmid + 1
    iflag = n - 2 * nmid
    allocate(xc(n), source = x)
    
    ! Sort the array in ascending order
    call sort(xc, .true.)

    ! Find the median
    if (iflag == 0) then
        rst = half * xc(nmid) + half * xc(nmidp1)
    else
        rst = xc(nmidp1)
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
    !! - <a href="https://en.wikipedia.org/wiki/Quantile" target="_blank">Wikipedia</a>
    real(real64), intent(in) :: x(:)
        !! An N-element array containing the data.
    real(real64), intent(in) :: q
        !! The quantile to compute (e.g. 0.25 computes the 25% quantile).
    real(real64) :: rst
        !! Quantile; empty/nonfinite data or q outside [0, 1] give NaN.

    ! Parameters
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    real(real64) :: a, b
    integer(int32) :: n, ib
    real(real64), allocatable :: xc(:)

    ! Initialization
    n = size(x)
    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (n == 0 .or. .not.all(ieee_is_finite(x))) return
    if (.not.ieee_is_finite(q) .or. q < 0.0d0 .or. q > 1.0d0) return
    allocate(xc(n), source = x)
    
    ! Sort the array in ascending order
    call sort(xc, .true.)

    ! Process
    a = (n + one) * q
    if (a <= one) then
        rst = xc(1)
        return
    end if
    if (a >= real(n, real64)) then
        rst = xc(n)
        return
    end if
    b = mod(a, one)
    ib = int(a - b, int32)
    
    ! Clamp index to valid range [1, n]
    ib = max(1, min(ib, n))
    
    if (ib >= n) then
        rst = xc(n)
    else
        rst = (one - b) * xc(ib) + b * xc(ib + 1)
    end if
    
    deallocate(xc)
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
        !! Trimmed mean; invalid p or empty/nonfinite data give NaN.

    ! Local Variables
    integer(int32) :: i1, i2, n
    real(real64) :: pv

    ! Initialization
    if (present(p)) then
        pv = p
    else
        pv = 0.05d0
    end if
    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (.not.ieee_is_finite(pv) .or. pv < 0.0d0 .or. pv >= 0.5d0) return
    if (size(x) == 0 .or. .not.all(ieee_is_finite(x))) return

    ! Sort the array into ascending order
    call sort(x, .true.)

    ! Find the limiting indices
    n = size(x)
    i1 = min(n, floor(n * pv, int32) + 1)
    i2 = max(1, n - floor(n * pv, int32))
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
        !! Sample covariance; fewer than two observations or nonfinite inputs
        !! give NaN. An unrepresentable result gives signed infinity.

    integer(int32) :: i, n
    real(real64) :: meanX, meanY, scaleX, scaleY, compensation, term, updated, log_magnitude
    real(real64), allocatable, dimension(:) :: shiftedX, shiftedY

    n = size(x)
    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (n < 2 .or. .not.all(ieee_is_finite(x)) .or. .not.all(ieee_is_finite(y))) return
    scaleX = maxval(abs(x))
    scaleY = maxval(abs(y))
    rst = 0.0d0
    if (scaleX == 0.0d0 .or. scaleY == 0.0d0) return
    shiftedX = shifted_scaled_values(x, scaleX)
    shiftedY = shifted_scaled_values(y, scaleY)
    meanX = mean(shiftedX)
    meanY = mean(shiftedY)
    compensation = 0.0d0
    do i = 1, n
        term = (shiftedX(i) - meanX) * (shiftedY(i) - meanY)
        updated = rst + term
        if (abs(rst) >= abs(term)) then
            compensation = compensation + ((rst - updated) + term)
        else
            compensation = compensation + ((term - updated) + rst)
        end if
        rst = updated
    end do
    rst = (rst + compensation) / real(n - 1, real64)
    if (rst == 0.0d0) return
    log_magnitude = log(abs(rst)) + log(scaleX) + log(scaleY)
    if (log_magnitude > log(huge(1.0d0))) then
        rst = sign(ieee_value(0.0d0, ieee_positive_inf), rst)
    else
        rst = sign(exp(log_magnitude), rst)
    end if
end function

! ------------------------------------------------------------------------------
pure function pooled_variance_1(si, ni) result(rst)
    !! Computes the pooled estimate of variance.
    real(real64), intent(in), dimension(:) :: si
        !! An N-element array containing the estimates for each of the N
        !! variances.
    integer(int32), intent(in), dimension(size(si)) :: ni
        !! An N-element array containing the number of data points in each
        !! of the data sets used to compute the variances in si.
    real(real64) :: rst
        !! Pooled variance; empty inputs, nonfinite/negative variances, or
        !! groups with fewer than two observations give NaN.

    integer(int32) :: i, k
    integer(int64) :: degrees_of_freedom
    real(real64) :: scale_value, term, compensation, updated

    k = size(si)
    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (k == 0 .or. any(ni < 2)) return
    if (.not.all(ieee_is_finite(si)) .or. any(si < 0.0d0)) return
    degrees_of_freedom = sum(int(ni, int64) - 1_int64)
    scale_value = maxval(si)
    rst = 0.0d0
    if (scale_value == 0.0d0) return
    compensation = 0.0d0
    do i = 1, k
        term = (si(i) / scale_value) * &
            (real(ni(i) - 1, real64) / real(degrees_of_freedom, real64)) - compensation
        updated = rst + term
        compensation = (updated - rst) - term
        rst = updated
    end do
    rst = scale_value * rst
end function

pure function pooled_variance_2(x) result(rst)
    !! Computes the pooled estimate of variance.
    type(array_container), intent(in), dimension(:) :: x
        !! An array of arrays of data.
    real(real64) :: rst
        !! Pooled variance; empty/unallocated groups or groups with fewer than
        !! two finite observations give NaN.

    integer(int32) :: i, k
    real(real64), allocatable, dimension(:) :: group_variances
    integer(int32), allocatable, dimension(:) :: group_sizes

    k = size(x)
    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (k == 0) return
    allocate(group_variances(k), group_sizes(k))
    do i = 1, k
        if (.not.allocated(x(i)%x)) return
        group_sizes(i) = size(x(i)%x)
        group_variances(i) = variance(x(i)%x)
    end do
    rst = pooled_variance_1(group_variances, group_sizes)
end function

! ------------------------------------------------------------------------------
end module
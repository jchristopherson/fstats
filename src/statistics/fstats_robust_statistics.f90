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

module fstats_robust_statistics
    use iso_fortran_env
    use fstats_descriptive_statistics
    implicit none
    private
    public :: median_absolute_deviation
    public :: tukey_biweight_psi
    public :: tukey_biweight_rho
    public :: tukey_biweight_irls_weight

contains
! ------------------------------------------------------------------------------
pure function median_absolute_deviation(x) result(rst)
    !! Computes the median absolute deviation of a data set by determining the
    !! median of distances of each data point from the data set median.
    real(real64), intent(in), dimension(:) :: x
        !! The data set.
    real(real64) :: rst
        !! The median absolute deviation of x.

    ! Local Variables
    integer(int32) :: n, i
    real(real64) :: x_median, x_abs, median_abs, max_value
    real(real64), allocatable, dimension(:) :: distances

    ! Process
    n = size(x)
    allocate(distances(n))
    x_median = median(x)
    max_value = huge(1.0d0)
    median_abs = abs(x_median)
    do i = 1, n
        if ((x(i) < 0.0d0 .and. x_median > 0.0d0) .or. &
            (x(i) > 0.0d0 .and. x_median < 0.0d0)) then
            x_abs = abs(x(i))
            if (x_abs > max_value - median_abs) then
                distances(i) = max_value
            else
                distances(i) = x_abs + median_abs
            end if
        else
            distances(i) = abs(x(i) - x_median)
        end if
    end do
    rst = median(distances)
end function

! ------------------------------------------------------------------------------
pure elemental function tukey_biweight_psi(u, k) result(rst)
    !! Computes Tukey's biweight score function for a standardized residual.
    !! The positive tuning constant k defaults to 4.685.
    !! $$
    !! \psi(u) =
    !! \begin{cases}
    !! u \left[1 - \left(\frac{u}{k}\right)^2\right]^2, & |u| \le k, \\
    !! 0, & |u| > k.
    !! \end{cases}
    !! $$
    real(real64), intent(in) :: u
        !! The standardized residual.
    real(real64), intent(in), optional :: k
        !! The positive cutoff/tuning constant; defaults to 4.685.
    real(real64) :: rst
        !! The score value psi(u).

    ! Local Variables
    real(real64) :: c
    
    ! Initialization
    if (present(k)) then
        c = k
    else
        c = 4.685d0
    end if

    ! Process
    if (abs(u) <= c) then
        rst = u * (1.0d0 - (u / c)**2)**2
    else
        rst = 0.0d0
    end if
end function

! ------------------------------------------------------------------------------
pure elemental function tukey_biweight_rho(u, k) result(rst)
    !! Computes Tukey's biweight rho loss for a standardized residual.
    !! The positive tuning constant k defaults to 4.685.
    !! $$
    !! \rho(u) =
    !! \begin{cases}
    !! \frac{k^2}{6} \left\{1 - \left[1 - \left(\frac{u}{k}\right)^2\right]^3\right\},
    !!     & |u| \le k, \\
    !! \frac{k^2}{6}, & |u| > k.
    !! \end{cases}
    !! $$
    real(real64), intent(in) :: u
        !! The standardized residual.
    real(real64), intent(in), optional :: k
        !! The positive cutoff/tuning constant; defaults to 4.685.
    real(real64) :: rst
        !! The loss value rho(u).

    ! Local Variables
    real(real64) :: c

    ! Initialization
    if (present(k)) then
        c = k
    else
        c = 4.685d0
    end if

    ! Process
    if (abs(u) <= c) then
        rst = (c**2 / 6.0d0) * (1.0d0 - (1.0d0 - (u / c)**2)**3)
    else
        rst = c**2 / 6.0d0
    end if
end function

! ------------------------------------------------------------------------------
pure elemental function tukey_biweight_irls_weight(u, k) result(rst)
    !! Computes Tukey's biweight IRLS weight for a standardized residual.
    !! The positive tuning constant k defaults to 4.685.
    !! $$
    !! w(u) =
    !! \begin{cases}
    !! \left[1 - \left(\frac{u}{k}\right)^2\right]^2, & |u| \le k, \\
    !! 0, & |u| > k.
    !! \end{cases}
    !! $$
    real(real64), intent(in) :: u
        !! The standardized residual.
    real(real64), intent(in), optional :: k
        !! The positive cutoff/tuning constant; defaults to 4.685.
    real(real64) :: rst
        !! The nonnegative IRLS weight.

    ! Local Variables
    real(real64) :: c

    ! Initialization
    if (present(k)) then
        c = k
    else
        c = 4.685d0
    end if

    ! Process
    if (abs(u) <= c) then
        rst = (1.0d0 - (u / c)**2)**2
    else
        rst = 0.0d0
    end if
end function

! ------------------------------------------------------------------------------
end module
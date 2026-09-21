! Copyright (c) 2023-2026 Jason Christopherson
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

program example
    use iso_fortran_env
    use fstats
    use fplot_core
    implicit none

    ! Description:
    ! This example uses rejection sampling to sample an F distribution.

    ! Local Variables
    integer(int32), parameter :: npts = 10000
    real(real64), parameter :: mu = 0.0d0
    real(real64), parameter :: sigma = 1.0d0
    real(real64), parameter :: xmin = 1.0d-4
    real(real64), parameter :: xmax = 2.0d0
    real(real64) :: x(npts)
    type(log_normal_distribution) :: dist

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_histogram) :: pd

    ! Perform the sampling
    dist%mean_value = mu
    dist%standard_deviation = sigma
    x = rejection_sample(dist, npts, xmin, xmax)

    ! Plot the resulting distribution
    call plt%initialize()
    call pd%define_data(x)
    call plt%push(pd)
    call plt%draw()
end program
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

module fstats_helper_routines
    use iso_fortran_env
    implicit none
    private
    public :: difference
    public :: factorial
contains

! ------------------------------------------------------------------------------
pure function difference(x) result(rst)
    !! Computes the difference between elements in an array.
    real(real64), intent(in), dimension(:) :: x
        !! The N-element array on which to operate.
    real(real64), allocatable, dimension(:) :: rst
        !! The (N-1)-element array containing the differences between adjacent
        !! elements.

    ! Local Variables
    integer(int32) :: i, n

    ! Process
    n = size(x)
    allocate(rst(n-1))
    do i = 1, n - 1
        rst(i) = x(i+1) - x(i)
    end do
end function

! ------------------------------------------------------------------------------
pure elemental function factorial(x) result(rst)
    !! Computes the factorial of X.
    real(real64), intent(in) :: x
        !! The value whose factorial is to be computed.
    real(real64) :: rst
        !! The result.
    rst = gamma(x + 1.0d0)
end function

! ------------------------------------------------------------------------------
end module
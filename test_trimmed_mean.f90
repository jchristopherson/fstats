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

program test_trimmed_mean
    use fstats
    use linalg, only : sort
    implicit none
    
    integer(int32), parameter :: n = 100
    real(real64) :: x(n), tm, ans
    integer(int32) :: i1, i2, i
    real(real64), parameter :: p = 0.05d0
    
    ! Create array with values 1 to 100
    do i = 1, n
        x(i) = real(i, real64)
    end do
    
    ! Manual calculation
    ! With p = 0.05, we remove 5% from bottom and 5% from top
    ! floor(100 * 0.05) = 5, so remove 5 from each end
    i1 = floor(real(n, real64) * p) + 1  ! Should be 6
    i2 = n - floor(real(n, real64) * p)  ! Should be 95
    
    print *, "n = ", n
    print *, "p = ", p
    print *, "i1 (correct) = ", i1
    print *, "i2 (correct) = ", i2
    print *, "Elements in range = ", i2 - i1 + 1
    
    ans = sum(x(i1:i2)) / (i2 - i1 + 1)
    print *, "Expected trimmed mean = ", ans
    
    ! Test function
    tm = trimmed_mean(x, p = p)
    print *, "Function returned = ", tm
    print *, "Match = ", abs(tm - ans) < 1.0d-6
    
end program test_trimmed_mean

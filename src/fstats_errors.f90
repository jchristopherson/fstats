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

! A module providing a set of routines to handle errors for the FSTATS library.
module fstats_errors
    use iso_fortran_env, only : int32
    implicit none

! ******************************************************************************
! ERROR CODES
! ------------------------------------------------------------------------------
    integer(int32), parameter :: FS_NO_ERROR = 0
    integer(int32), parameter :: FS_ARRAY_SIZE_ERROR = 10000
    integer(int32), parameter :: FS_MATRIX_SIZE_ERROR = 10001
    integer(int32), parameter :: FS_INVALID_INPUT_ERROR = 10002
    integer(int32), parameter :: FS_MEMORY_ERROR = 10003
    integer(int32), parameter :: FS_UNDERDEFINED_PROBLEM_ERROR = 10004
    integer(int32), parameter :: FS_TOLERANCE_TOO_SMALL_ERROR = 10005
    integer(int32), parameter :: FS_TOO_FEW_ITERATION_ERROR = 10006
    integer(int32), parameter :: FS_INVALID_ARGUMENT_ERROR = 10007
    integer(int32), parameter :: FS_ZERO_VARIANCE_ERROR = 10008
    integer(int32), parameter :: FS_NULL_POINTER_ERROR = 10009
    integer(int32), parameter :: FS_POLYNOMIAL_ORDER_ERROR = 10010
    integer(int32), parameter :: FS_NONMONOTONIC_ARRAY_ERROR = 10011
    integer(int32), parameter :: FS_UNINITIALIZED_OBJECT_ERROR = 10012

end module
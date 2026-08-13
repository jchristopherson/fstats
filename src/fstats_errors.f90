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
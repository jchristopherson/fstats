! A module providing a set of routines to handle errors for the FSTATS library.
module fstats_errors
    use ferror
    use iso_fortran_env, only : int32
    implicit none
    private
    public :: FS_NO_ERROR
    public :: FS_ARRAY_SIZE_ERROR
    public :: FS_INVALID_INPUT_ERROR
    public :: FS_OUT_OF_MEMORY_ERROR
    public :: FS_UNDERDEFINED_PROBLEM_ERROR
    public :: FS_TOLERANCE_TOO_SMALL_ERROR

! ******************************************************************************
! ERROR CODES
! ------------------------------------------------------------------------------
    integer(int32), parameter :: FS_NO_ERROR = 0
    integer(int32), parameter :: FS_ARRAY_SIZE_ERROR = 10000
    integer(int32), parameter :: FS_INVALID_INPUT_ERROR = 10001
    integer(int32), parameter :: FS_MEMORY_ERROR = 10002
    integer(int32), parameter :: FS_UNDERDEFINED_PROBLEM_ERROR = 10003
    integer(int32), parameter :: FS_TOLERANCE_TOO_SMALL_ERROR = 10004
    integer(int32), parameter :: FS_TOO_FEW_ITERATION_ERROR = 10005

! ------------------------------------------------------------------------------
    integer(int32), parameter :: MESSAGE_SIZE = 1024

contains
! ------------------------------------------------------------------------------
    !> @brief Reports a memory allocation related error.
    !!
    !! @param[in,out] err The error handling object.
    !! @param[in] fname The name of the routine in which the error occurred.
    !! @param[in] code The error code returned by the allocation routine.
    subroutine report_memory_error(err, fname, code)
        ! Arguments
        class(errors), intent(inout) :: err
        character(len = *), intent(in) :: fname
        integer(int32), intent(in) :: code

        ! Variables
        character(len = MESSAGE_SIZE) :: msg

        ! Process
        write(msg, 100) &
            "A memory allocation error occurred with code ", code, "."
        call err%report_error(fname, trim(msg), FS_MEMORY_ERROR)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Reports an array size error.
    !!
    !! @param[in,out] err The error handling object.
    !! @param[in] fname The name of the routine in which the error occurred.
    !! @param[in] name The name of the array.
    !! @param[in] expect The expected size of the array.
    !! @param[in] actual The actual size of the array.
    subroutine report_array_size_error(err, fname, name, expect, actual)
        ! Arguments
        class(errors), intent(inout) :: err
        character(len = *), intent(in) :: fname, name
        integer(int32), intent(in) :: expect, actual

        ! Variables
        character(len = MESSAGE_SIZE) :: msg

        ! Process
        write(msg, 100) "Expected array " // name // " to be of length ", &
            expect, ", but found it to be of length ", actual, "."
        call err%report_error(fname, trim(msg), FS_ARRAY_SIZE_ERROR)

        ! Formatting
100     format(A, I0, A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Reports an error relating to two arrays not being the same size
    !! when they should be the same size.
    !!
    !! @param[in,out] The error handling object.
    !! @param[in] fname The name of the routine in which the error occurred.
    !! @param[in] name1 The name of the first array.
    !! @param[in] name2 The name of the second array.
    !! @param[in] size1 The size of the first array.
    !! @param[in] size2 The size of the second array.
    subroutine report_arrays_not_same_size_error(err, fname, name1, name2, &
        size1, size2)
        ! Arguments
        class(errors), intent(inout) :: err
        character(len = *), intent(in) :: fname, name1, name2
        integer(int32), intent(in) :: size1, size2

        ! Local Variables
        character(len = MESSAGE_SIZE) :: msg

        ! Process
        write(msg, 100) "Array " // name1 // " and array " // name2 // &
            "were expected to be the same size, but instead were found " // &
            "to be sized ", size1, " and ", size2 " respectively."
        call err%report_error(fname, trim(msg), FS_ARRAY_SIZE_ERROR)

        ! Formatting
100     format(A, I0, A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
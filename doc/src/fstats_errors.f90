! A module providing a set of routines to handle errors for the FSTATS library.
module fstats_errors
    use ferror
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

! ------------------------------------------------------------------------------
    integer(int32), private, parameter :: MESSAGE_SIZE = 1024

contains
! ------------------------------------------------------------------------------
    subroutine report_memory_error(err, fname, code)
        !! Reports a memory allocation related error.
        class(errors), intent(inout) :: err
            !! The error handling object.
        character(len = *), intent(in) :: fname
            !! The name of the routine in which the error occurred.
        integer(int32), intent(in) :: code
            !! The error code returned by the allocation routine.

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
    subroutine report_array_size_error(err, fname, name, expect, actual)
        !! Reports an array size error.
        class(errors), intent(inout) :: err
            !! The error handling object.
        character(len = *), intent(in) :: fname
            !! The name of the routine in which the error occurred.
        character(len = *), intent(in) :: name
            !! The name of the array.
        integer(int32), intent(in) :: expect
            !! The expected size of the array.
        integer(int32), intent(in) :: actual
            !! The actual size of the array.

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
    subroutine report_matrix_size_error(err, fname, name, expect_rows, &
        expect_cols, actual_rows, actual_cols)
        !! Reports a matrix size error.
        class(errors), intent(inout) :: err
            !! The error handling object.
        character(len = *), intent(in) :: fname
            !! The name of the routine in which the error occurred.
        character(len = *), intent(in) :: name
            !! The name of the matrix.
        integer(int32), intent(in) :: expect_rows
            !! The expected number of rows.
        integer(int32), intent(in) :: expect_cols
            !! The expected number of columns.
        integer(int32), intent(in) :: actual_rows
            !! The actual number of rows.
        integer(int32), intent(in) :: actual_cols
            !! The actual number of columns.
        
        ! Variables
        character(len = MESSAGE_SIZE) :: msg

        ! Process
        write(msg, 100) "Expected matrix " // name // " to be of size (", &
            expect_rows, ", ", expect_cols, "), but found it to be of size (", &
            actual_rows, ", ", actual_cols, ")."
        call err%report_error(fname, trim(msg), FS_MATRIX_SIZE_ERROR)

        ! Formatting
100     format(A, I0, A, I0, A, I0, A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine report_arrays_not_same_size_error(err, fname, name1, name2, &
        size1, size2)
        !! Reports an error relating to two arrays not being the same size
        !! when they should be the same size.
        class(errors), intent(inout) :: err
            !! The error handling object.
        character(len = *), intent(in) :: fname
            !! The name of the routine in which the error occurred.
        character(len = *), intent(in) :: name1
            !! The name of the first array.
        character(len = *), intent(in) :: name2
            !! The name of the second array.
        integer(int32), intent(in) :: size1
            !! The size of the first array.
        integer(int32), intent(in) :: size2
            !! The size of the second array.

        ! Local Variables
        character(len = MESSAGE_SIZE) :: msg

        ! Process
        write(msg, 100) "Array " // name1 // " and array " // name2 // &
            "were expected to be the same size, but instead were found " // &
            "to be sized ", size1, " and ", size2, " respectively."
        call err%report_error(fname, trim(msg), FS_ARRAY_SIZE_ERROR)

        ! Formatting
100     format(A, I0, A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine report_underdefined_error(err, fname, expect, actual)
        !! Reports an underdefined problem error.
        class(errors), intent(inout) :: err
            !! The error handling object.
        character(len = *), intent(in) :: fname
            !! The name of the routine in which the error occurred.
        integer(int32), intent(in) :: expect
            !! The expected minimum number of equations.
        integer(int32), intent(in) :: actual
            !! The actual number of equations.

        ! Local Variables
        character(len = MESSAGE_SIZE) :: msg

        ! Process
        write(msg, 100) "The problem is underdefined.  The number of " // &
            "equations was found to be ", actual, &
            ", but must be at least equal to the number of unknowns ", &
            expect, "."
        call err%report_error(fname, trim(msg), FS_UNDERDEFINED_PROBLEM_ERROR)
        
        ! Formatting
100     format(A, I0, A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine report_iteration_count_error(err, fname, msg, mincount)
        !! Reports an iteration count error.
        class(errors), intent(inout) :: err
            !! The error handling object.
        character(len = *), intent(in) :: fname
            !! The name of the routine in which the error occurred.
        character(len = *), intent(in) :: msg
            !! The error message.
        integer(int32), intent(in) :: mincount
            !! The minimum iteration count expected.

        ! Local Variables
        character(len = MESSAGE_SIZE) :: emsg

        ! Process
        write(emsg, 100) msg // "  A minimum of ", mincount, " is expected."
        call err%report_error(fname, trim(emsg), FS_TOO_FEW_ITERATION_ERROR)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine report_zero_variance_warning(err, fname)
        !! Issues a warning when a zero variance is encountered in an unexpected
        !! place.
        class(errors), intent(inout) :: err
            !! The error handling object.
        character(len = *), intent(in) :: fname
            !! The name of the routine in which the error occurred.

        ! Process
        call err%report_warning(fname, &
            "A zero-valued variance was encountered when not expected.", &
            FS_ZERO_VARIANCE_ERROR)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine report_null_pointer_error(err, fname, name)
        !! Reports a null pointer error.
        class(errors), intent(inout) :: err
            !! The error handling object.
        character(len = *), intent(in) :: fname
            !! The name of the routine in which the error occurred.
        character(len = *), intent(in) :: name
            !! The parameter name.

        ! Process
        call err%report_error(fname, &
            "The parameter " // name // " was found to be null.", &
            FS_NULL_POINTER_ERROR)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine report_polynomial_order_error(err, fname, order, minorder)
        !! Reports an error related to the order of the polynomial.
        class(errors), intent(inout) :: err
            !! The error handling object.
        character(len = *), intent(in) :: fname
            !! The name of the routine in which the error occurred.
        integer(int32), intent(in) :: order
            !! The supplied order.
        integer(int32), intent(in) :: minorder
            !! The minimum required polynomial order.

        ! Local Variables
        character(len = 256) :: msg

        ! Process
        write(msg, 100) "The supplied polynomial order of ", order, &
            " is below the minimum required order of ", minorder, "."
        call err%report_error(fname, trim(msg), FS_POLYNOMIAL_ORDER_ERROR)

        ! Formatting
100     format(A, I0, A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine report_nonmonotonic_array_error(err, fname, name)
        !! Reports an error related to a nonmonotonic array when a monotonic
        !! array was expected.
        class(errors), intent(inout) :: err
            !! The error handling object.
        character(len = *), intent(in) :: fname
            !! The name of the routine in which the error occurred.
        character(len = *), intent(in) :: name
            !! The parameter name.

        ! Local Variables
        character(len = 256) :: msg

        ! Process
        write(msg, 100) "Array ", name, &
            " was expected to be monotonic, but was found to be nonmonotonic."
        call err%report_error(fname, trim(msg), FS_NONMONOTONIC_ARRAY_ERROR)

        ! Formatting
100     format(A, A, A)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine report_uninitialized_object_error(err, fname, obj)
        !! Reports an uninitialized object error.
        class(errors), intent(inout) :: err
            !! The error handling object.
        character(len = *), intent(in) :: fname
            !! The name of the routine in which the error occurred.
        character(len = *), intent(in) :: obj
            !! The name of the object or type.

        ! Local Variables
        character(len = 256) :: msg

        ! Process
        write(msg, 100) "Object ", obj, " is uninitialized."
        call err%report_error(fname, trim(msg), FS_UNINITIALIZED_OBJECT_ERROR)

        ! Formatting
100     format(A, A, A)
    end subroutine

! ------------------------------------------------------------------------------
end module
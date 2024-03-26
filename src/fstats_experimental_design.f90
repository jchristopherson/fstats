module fstats_experimental_design
    use iso_fortran_env
    use fstats_errors
    implicit none
    private
    public :: get_full_factorial_matrix_size
    public :: full_factorial
contains
! ------------------------------------------------------------------------------
subroutine get_full_factorial_matrix_size(vars, m, n, err)
    !! Computes the appropriate size for a full-factorial design table.
    integer(int32), intent(in) :: vars(:)
        !! An M-element array containing the M factors to study.  Each 
        !! of the M entries to the array is expected to contain the 
        !! number of options for that particular factor to explore.  
        !! This value must be greater than or equal to 1.
    integer(int32), intent(out) :: m
        !! The number of rows for the table.
    integer(int32), intent(out) :: n
        !! The number of columns for the table.
    class(errors), intent(inout), optional, target :: err
        !! A mechanism for communicating errors and warnings to the 
        !! caller.  Possible warning and error codes are as follows.
        !! - FS_NO_ERROR: No errors encountered.
        !! - FS_INVALID_INPUT_ERROR: Occurs if any items in vars are 
        !!      less than 1.

    ! Local Variables
    integer(int32) :: i
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = 256) :: errmsg
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    m = 0
    n = 0

    ! Ensure every value is greater than 1
    do i = 1, size(vars)
        if (vars(i) < 1) then
            write(errmsg, 100) "A value less than 1 was found at index ", &
                i, " of the input array.  All values must be greater " // &
                "than or equal to 1."
            call errmgr%report_error("get_full_factorial_matrix_size", &
                trim(errmsg), FS_INVALID_INPUT_ERROR)
            return
        end if
    end do

    ! Process
    m = product(vars)
    n = size(vars)

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
subroutine full_factorial(vars, tbl, err)
    !! Computes a table with values scaled from 1 to N describing a 
    !! full-factorial design.
    !!
    !! ```fortran
    !! program example
    !!     use iso_fortran_env
    !!     use fstats
    !!     implicit none
    !!
    !!     ! Local Variables
    !!     integer(int32) :: i, vars(3), tbl(24, 3)
    !!
    !!     ! Define the number of design points for each of the 3 factors to study
    !!     vars = [2, 4, 3]
    !!
    !!     ! Determine the design table
    !!     call full_factorial(vars, tbl)
    !!
    !!     ! Display the table
    !!     do i = 1, size(tbl, 1)
    !!         print *, tbl(i,:)
    !!     end do
    !! end program
    !! ```
    !! The above program produces the following output.
    !! ```text
    !! 1           1           1
    !! 1           1           2
    !! 1           1           3
    !! 1           2           1
    !! 1           2           2
    !! 1           2           3
    !! 1           3           1
    !! 1           3           2
    !! 1           3           3
    !! 1           4           1
    !! 1           4           2
    !! 1           4           3
    !! 2           1           1
    !! 2           1           2
    !! 2           1           3
    !! 2           2           1
    !! 2           2           2
    !! 2           2           3
    !! 2           3           1
    !! 2           3           2
    !! 2           3           3
    !! 2           4           1
    !! 2           4           2
    !! 2           4           3
    !! ```
    integer(int32), intent(in) :: vars(:)
        !! An M-element array containing the M factors to study.  
        !! Each of the M entries to the array is expected to contain 
        !! the number of options for that particular factor to explore. 
        !! This value must be greater than or equal to 1.
    integer(int32), intent(out) :: tbl(:,:)
        !! A table where the design will be written.  Use 
        !! get_full_factorial_matrix_size to determine the appropriate 
        !! table size.
    class(errors), intent(inout), optional, target :: err
        !! A mechanism for communicating errors and warnings to the 
        !! caller.  Possible warning and error codes are as follows.
        !! - FS_NO_ERROR: No errors encountered.
        !! - FS_INVALID_INPUT_ERROR: Occurs if any items in vars are 
        !!      less than 1.
        !! - FS_ARRAY_SIZE_ERROR: Occurs if tbl is not properly sized.

    ! Local Variables
    integer(int32) :: i, col, stride, last, val, m, n
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Verify the size of the input table
    call get_full_factorial_matrix_size(vars, m, n, errmgr)
    if (errmgr%has_error_occurred()) return
    if (size(tbl, 1) /= m .or. size(tbl, 2) /= n) then
        call report_matrix_size_error(errmgr, "full_factorial", &
            "tbl", m, n, size(tbl, 1), size(tbl, 2))
        return
    end if

    ! Process
    do col = 1, n
        stride = 1
        if (col /= n) stride = product(vars(col+1:n))
        val = 1
        do i = 1, m, stride
            last = i + stride - 1
            tbl(i:last,col) = val
            val = val + 1
            if (val > vars(col)) val = 1
        end do
    end do
end subroutine

! ------------------------------------------------------------------------------
end module
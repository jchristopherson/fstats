submodule (fstats) experimental_design_implementation
    use fstats_errors
    implicit none
contains
! ------------------------------------------------------------------------------
module subroutine get_full_factorial_matrix_size_int32(vars, m, n, err)
    ! Arguments
    integer(int32), intent(in) :: vars(:)
    integer(int32), intent(out) :: m, n
    class(errors), intent(inout), optional, target :: err

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
            call errmgr%report_error("get_full_factorial_matrix_size_int32", &
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
module subroutine full_factorial_int32(vars, tbl, err)
    ! Arguments
    integer(int32), intent(in) :: vars(:)
    integer(int32), intent(out) :: tbl(:,:)
    class(errors), intent(inout), optional, target :: err

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
        call report_matrix_size_error(errmgr, "full_factorial_int32", &
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
end submodule
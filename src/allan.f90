submodule (fstats) allan
    use fstats_errors
contains
! ------------------------------------------------------------------------------
module function allan_variance(x, err) result(rst)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable, dimension(:) :: rst

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    integer(int32) :: flag, j, m, n, limit, nr
    real(real64), allocatable, dimension(:) :: tall1, tall2
    real(real64) :: temp
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Initialization
    n = size(x)
    limit = n
    nr = floor(0.5 * n) - 1
    allocate(tall1(n - 1), source = x(:n-1), stat = flag)
    if (flag == 0) allocate(tall2(n - 1), source = x(2:n))
    if (flag == 0) allocate(rst(nr), source = 0.0d0)
    if (flag /= 0) go to 10

    ! Process
    do m = 1, nr
        temp = sum((tall2 - tall1)**2)
        do j = 1, limit - 1
            tall1(j) = tall1(j) + x(m + j)
            tall2(j) = tall2(j + 1) + x(2 * m + j + 1)
        end do
        limit = limit - 2
        rst(m) = temp / (2.0d0 * (n - 2 * m + 1) * m**2)
    end do


    ! End
    return

    ! Memory Error Handling
10  continue
    call report_error(errmgr, "allan_variance", flag)
    return
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end submodule
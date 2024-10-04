module fstats_experimental_design
    use iso_fortran_env
    use fstats_errors
    implicit none
    private
    public :: get_full_factorial_matrix_size
    public :: full_factorial
    public :: doe_fit_model
    public :: doe_evaluate_model
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
subroutine doe_fit_model(nway, x, y, beta, err)
    !! Fits a Taylor series model to the provided data.
    !!
    !! $$ Y = \beta_{0} + \sum_{i=1}^{n} \beta_{i} X_{i} + \sum_{i=1}^{n} 
    !! \sum_{j=1 \\ i \neq j}^{n} \beta_{ij} X_{i} X_{j} + \sum_{i=1}^{n} 
    !! \sum_{j=1}^{n} \sum_{k=1 \\ i \neq j \neq k}^{n} \beta_{ijk} X_{i} 
    !! X_{j} X_{k} + ... $$
    integer(int32), intent(in) :: nway
        !! The number of interaction levels.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! used to produce the results.
    real(real64), intent(in), dimension(:) :: y
        !! An M-element array containing the results from the M experiments.
    real(real64), intent(out), allocatable, dimension(:) :: beta
        !! An allocatable array that will contain the fitted model parameters.
    class(errors), intent(inout), optional, target :: err
        !! A mechanism for communicating errors and warnings to the 
        !! caller.  Possible warning and error codes are as follows.
        !! - FS_NO_ERROR: No errors encountered.
        !! - FS_ARRAY_SIZE_ERROR: Occurs if x and y are not properly sized
        !!      relative to one another.
end subroutine

! ------------------------------------------------------------------------------
function doe_evaluate_model(nway, beta, x, map, err) result(rst)
    !! Evaluates the model of the following form.
    !!
    !! $$ Y = \beta_{0} + \sum_{i=1}^{n} \beta_{i} X_{i} + \sum_{i=1}^{n} 
    !! \sum_{j=1 \\ i \neq j}^{n} \beta_{ij} X_{i} X_{j} + \sum_{i=1}^{n} 
    !! \sum_{j=1}^{n} \sum_{k=1 \\ i \neq j \neq k}^{n} \beta_{ijk} X_{i} 
    !! X_{j} X_{k} + ... $$
    integer(int32), intent(in) :: nway
        !! The number of interaction levels.  Currently, this algorithm supports
        !! a maximum of three-way interaction.
    real(real64), intent(in), dimension(:) :: beta
        !! The model coefficients.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! at which to evaluate the model.
    logical, intent(in), optional, target, dimension(:) :: map
        !! An optional array of the same size as beta that can be used to
        !! eliminate a parameter from the model (false), or keep a parameter
        !! in the model (true).  If not supplied, all parameters will be assumed
        !! to be part of the model as if the array were filled with all true
        !! values.
    class(errors), intent(inout), optional, target :: err
        !! A mechanism for communicating errors and warnings to the 
        !! caller.  Possible warning and error codes are as follows.
        !! - FS_NO_ERROR: No errors encountered.
        !! - FS_ARRAY_SIZE_ERROR: Occurs if beta and map are not properly sized
        !!      relative to one another.
    real(real64), allocatable, dimension(:) :: rst
        !! The resulting M-element array.

    ! Local Variables
    integer(int32) :: m, flag
    logical, pointer, dimension(:) :: mapptr
    logical, allocatable, target, dimension(:) :: nmap

    ! Initialization
    m = size(x, 1)

    ! Input Checking
    if (nway < 1 .or. nway > 3) then
        ! TO DO: Error - must be at least 1, but not more than 3
    end if
    ! TO DO: Check the size of beta

    ! Memory Allocations
    allocate(rst(m), stat = flag)
    if (flag /= 0) then
        ! TO DO: Error - memory issue
    end if

    ! Set up the map parameters
    if (present(map)) then
        ! TO DO: Check the size of map
        mapptr => map
    else
        allocate(nmap(size(beta)), stat = flag, source = .true.)
        if (flag /= 0) then
            ! TO DO: Error - memory issue
        end if
        mapptr => nmap
    end if

    ! Process
    call doe_eval_engine(nway, beta, x, mapptr, rst)
end function

! ----------
subroutine doe_eval_engine(nway, beta, x, map, y)
    ! Driver routine for "doe_evaluate_model" that performs the actual 
    ! calculations but forgoes any error checking.  This should not be exposed 
    ! as part of the public API.
    integer(int32), intent(in) :: nway
    real(real64), intent(in), dimension(:) :: beta
    real(real64), intent(in), dimension(:,:) :: x
    logical, intent(in), dimension(:) :: map
    real(real64), intent(out), dimension(:) :: y

    ! Local Variables
    integer(int32) :: i1, i2, n

    ! Initialization
    n = size(x, 2)
    if (map(1)) then
        y = beta(1)
    else
        y = 0.0d0
    end if

    ! Process
    if (nway >= 1) then
        i1 = 2
        i2 = i1 + n - 1
        call doe_eval_1(beta(i1:i2), x, map(i1:i2), y)
    end if
    if (nway >= 2) then
        i1 = i2 + 1
        i2 = i1 + n * (n - 1) - 1
        call doe_eval_2(beta(i1:i2), x, map(i1:i2), y)
    end if
    if (nway >= 3) then
        i1 = i2 + 1
        i2 = i1 + n * (n**2 - 1) - 1
        call doe_eval_3(beta(i1:i2), x, map(i1:i2), y)
    end if
end subroutine

! ----------
subroutine doe_eval_1(beta, x, map, y)
    !! Evaluates the main effect term.
    !!
    !! $$ Y = Y + /sum_{i=1}^{n} \beta_{i} X_{i} $$
    real(real64), intent(in), dimension(:) :: beta
        !! The model coefficients for just this portion of the model.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! at which to evaluate the model.
    logical, intent(in), dimension(:) :: map
        !! The usage map corresponding to the model coefficients for just this
        !! portion of the model.
    real(real64), intent(inout), dimension(:) :: y
        !! On input, an M-element array containing the existing portion of the 
        !! model.  On output, this array is updated to include the main effects.

    ! Local Variables
    integer(int32) :: i, n

    ! Initialization
    n = size(x, 2)

    ! Process
    do i = 1, n
        if (.not.map(i)) cycle
        y = y + beta(i) * x(:,i)
    end do
end subroutine

! ----------
subroutine doe_eval_2(beta, x, map, y)
    !! Evaluates the two-way interaction term.
    !!
    !! $$ Y = Y + /sum_{i=1}^{n} /sum_{j=1 // i /neq j}^{n} \beta_{i} X_{i} 
    !! X_{j} $$
    real(real64), intent(in), dimension(:) :: beta
        !! The model coefficients for just this portion of the model.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! at which to evaluate the model.
    logical, intent(in), dimension(:) :: map
        !! The usage map corresponding to the model coefficients for just this
        !! portion of the model.
    real(real64), intent(inout), dimension(:) :: y
        !! On input, an M-element array containing the existing portion of the 
        !! model.  On output, this array is updated to include the two-way
        !! interactions.

    ! Local Variables
    integer(int32) :: i, j, k, n

    ! Initialization
    n = size(x, 2)

    ! Process
    k = 1
    do i = 1, n
        do j = 1, n
            if (i == j) cycle
            if (.not.map(k)) cycle
            y = y + beta(k) * x(:,i) * x(:,j)
            k = k + 1
        end do
    end do
end subroutine

! ----------
subroutine doe_eval_3(beta, x, map, y)
    !! Evaluates the three-way interaction term.
    !!
    !! $$ Y = Y + /sum_{i=1}^{n} /sum_{j=1 // i /neq j}^{n} \beta_{i} X_{i} 
    !! X_{j} $$
    real(real64), intent(in), dimension(:) :: beta
        !! The model coefficients for just this portion of the model.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! at which to evaluate the model.
    logical, intent(in), dimension(:) :: map
        !! The usage map corresponding to the model coefficients for just this
        !! portion of the model.
    real(real64), intent(inout), dimension(:) :: y
        !! On input, an M-element array containing the existing portion of the 
        !! model.  On output, this array is updated to include the three-way
        !! interactions.

    ! Local Variables
    integer(int32) :: i, j, k, ii, n

    ! Initialization
    n = size(x, 2)

    ! Process
    ii = 1
    do i = 1, n
        do j = 1, n
            do k = 1, n
                if (i == j .and. j == k) cycle
                if (.not.map(ii)) cycle
                y = y + beta(ii) * x(:,i) * x(:,j) * x(:,k)
                ii = ii + 1
            end do
        end do
    end do
end subroutine

! ------------------------------------------------------------------------------
end module
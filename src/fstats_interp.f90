module fstats_interp
    use iso_fortran_env
    use fstats_errors
    use ferror
    implicit none
    private
    public :: linear_interpolator

    type base_interpolator
        !! A base object containing the functionallity to support 1D
        !! interpolation operations, regardless of the type of 1D interpolation.
        real(real64), public, allocatable, dimension(:) :: x
            !! An N-element array containing the raw x-coordinate data.
        real(real64), public, allocatable, dimension(:) :: y
            !! An N-element array containing the raw y-coordinate data.
        integer(int32), private :: m_order
            !! The order of the interpolating polynomial.
        integer(int32), private :: m_nearestIndex
            !! The saved nearest index from the last search.  Utilizing this
            !! value improves the behavior of the search.
        integer(int32), private :: m_method
            !! An integer used to identify the search method to use during
            !! interpolation.
        integer(int32), private :: m_indexDelta
            !! The bracket index difference.
    contains
        procedure, public :: initialize => bi_init
        procedure, public :: locate => bi_locate
        procedure, public :: hunt => bi_hunt
        procedure, public :: size => bi_size
    end type

! ------------------------------------------------------------------------------
    type :: linear_interpolator
        type(base_interpolator), private :: m_manager
            !! An object to manage the interpolation process.
    contains
        procedure, public :: initialize => li_init
        procedure, public :: interpolate => li_interp
        procedure, private :: raw_interp => li_raw_interp
    end type

contains
! ******************************************************************************
! BASE_INTERPOLATOR
! ------------------------------------------------------------------------------
subroutine bi_init(this, x, y, order, err)
    !! Initializes the interpolation object.
    class(base_interpolator), intent(inout) :: this
        !! The base_interpolator object.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the x-coordinate data in either 
        !! monotonically increasing or decreasing order.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the y-coordinate data.
    integer(int32), intent(in) :: order
        !! The order of the interpolating polynomial.
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.

    ! Local Variables
    integer(int32) :: i, n, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)

    ! Input Checking
    if (size(y) /= n) then
        call report_array_size_error(errmgr, "bi_init", "y", n, size(y))
        return
    end if

    if (order < 1) then
        ! TO DO: Report the error
        return
    end if

    ! TO DO: Ensure X is monotonic

    ! Process
    this%m_order = order
    this%m_nearestIndex = 1
    this%m_indexDelta = 1
    this%m_method = 0

    if (allocated(this%x)) deallocate(this%x)
    if (allocated(this%y)) deallocate(this%y)
    allocate(this%x(n), this%y(n), stat = flag)
    if (flag /= 0) go to 10
    do i = 1, n
        this%x(i) = x(i)
        this%y(i) = y(i)
    end do
    
    ! End
    return

    ! Memory Error Handling
10  continue
    call report_memory_error(errmgr, "bi_init", flag)
    return
end subroutine

! ------------------------------------------------------------------------------
function bi_locate(this, x) result(rst)
    !! Locates the index below the nearest entry to x.
    class(base_interpolator), intent(inout) :: this
        !! The base_interpolator object.
    real(real64), intent(in) :: x
        !! The x-coordinate of interest.
    integer(int32) :: rst
        !! The index below the nearest entry to x.  This value is always within
        !! the range of possible indices for the stored data.

    ! Local Variables
    logical :: ascnd
    integer(int32) :: ju, jm, jl, mm, n

    ! Initialization
    n = size(this%x)
    mm = this%m_order + 1
    ascnd = this%x(n) >= this%x(1)
    jl = 1
    ju = n

    ! Process
    do while (ju - jl > 1)
        jm = (ju + jl) / 2
        if ((x >= this%x(jm)) .eqv. ascnd) then
            jl = jm
        else
            ju = jm
        end if
    end do
    if (abs(jl - this%m_nearestIndex) > this%m_indexDelta) then
        this%m_method = 0
    else
        this%m_method = 1
    end if
    this%m_nearestIndex = jl

    ! Output
    if (x == this%x(1)) then
        rst = 1
    else if (x == this%x(n)) then
        rst = n - 1
    else
        rst = jl
    end if
end function

! ------------------------------------------------------------------------------
function bi_hunt(this, x) result(rst)
    !! A search routine, similar to locate, but tailored towards searching
    !! when the initial location estimate is very near the desired point.  If
    !! not near, this method is rather inefficient; however, when near, this
    !! method is quite efficient when compared with locate.
    class(base_interpolator), intent(inout) :: this
        !! The base_interpolator object.
    real(real64), intent(in) :: x
        !! The x-coordinate of interest.
    integer(int32) :: rst
        !! The index below the nearest entry to x.  This value is always within
        !! the range of possible indices for the stored data.

    ! Local Variables
    integer(int32) :: jl, jm, ju, inc, n, mm
    logical :: ascnd

    ! Initialization
    n = size(this%x)
    mm = this%m_order + 1
    jl = this%m_nearestIndex
    inc = 1
    ascnd = this%x(n) >= this%x(1)

    ! Process
    if (jl < 1 .or. jl > n) then
        jl = 1
        ju = n
    else
        if (x >= this%x(jl) .eqv. ascnd) then
            do
                ju = jl + inc
                if (ju >= n) then
                    ju = n
                    exit
                else if (x < this%x(ju) .eqv. ascnd) then
                    exit
                else
                    jl = ju
                    inc = inc + inc
                end if
            end do
        else
            ju = jl
            do
                jl = jl - inc
                if (jl <= 1) then
                    jl = 1
                    exit
                else if (x >= this%x(jl) .eqv. ascnd) then
                    exit
                else
                    ju = jl
                    inc = inc + inc
                end if
            end do
        end if
    end if

    ! THe hunt is done, so begin the final bisection
    do while (ju - jl > 1)
        jm = (ju + jl) / 2
        if (x >= this%x(jm) .eqv. ascnd) then
            jl = jm
        else
            ju = jm
        end if
    end do

    ! Check to see if we should hunt or locate next time around
    if (abs(jl - this%m_nearestIndex) > this%m_indexDelta) then
        this%m_method = 0
    else
        this%m_method = 1
    end if
    this%m_nearestIndex = jl

    ! Output
    if (x == this%x(1)) then
        rst = 1
    else if (x == this%x(n)) then
        rst = n - 1
    else
        rst = jl
    end if
end function

! ------------------------------------------------------------------------------
pure function bi_size(this) result(rst)
    !! Gets the size of the stored raw data array.
    class(base_interpolator), intent(in) :: this
        !! The base_interpolator object.
    integer(int32) :: rst
        !! The size of the stored data.

    if (allocated(this%x)) then
        rst = size(this%x)
    else
        rst = 0
    end if
end function

! ******************************************************************************
! LINEAR_INTERPOLATOR
! ------------------------------------------------------------------------------
subroutine li_init(this, x, y, err)
    !! Initializes the interpolation object.
    class(linear_interpolator), intent(inout) :: this
        !! The linear_interpolator object.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the x-coordinate data in either 
        !! monotonically increasing or decreasing order.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the y-coordinate data.
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Initialize the managing object
    call this%m_manager%initialize(x, y, 1, err = errmgr)
    if (errmgr%has_error_occurred()) return
end subroutine

! ------------------------------------------------------------------------------
subroutine li_interp(this, x, yi, err)
    !! Performs the interpolation.
    class(linear_interpolator), intent(inout) :: this
        !! The linear_interpolator object.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the x values at which to compute the
        !! interpolation.
    real(real64), intent(out), dimension(:) :: yi
        !! An N-element array  containing the interpolated data.
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.

    ! Local Variables
    integer(int32) :: i, j, n
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)

    ! Input Checking
    if (size(yi) /= n) then
        call report_array_size_error(errmgr, "li_interp", "yi", n, size(yi))
        return
    end if

    ! Perform the interpolation
    do i = 1, n
        if (this%m_manager%m_method == 1) then
            j = this%m_manager%hunt(x(i))
        else
            j = this%m_manager%locate(x(i))
        end if
        yi(j) = this%raw_interp(j, x(i))
    end do
end subroutine

! ------------------------------------------------------------------------------
function li_raw_interp(this, j, x) result(rst)
    !! Raw, single-value interpolation routine.
    class(linear_interpolator), intent(inout) :: this
        !! The linear_interpolator object.
    integer(int32), intent(in) :: j
        !! The array index immediately below the x array location needed.
    real(real64), intent(in) :: x
        !! The value at which to compute the interpolation.
    real(real64) :: rst
        !! The interpolated value.

    ! Process
    if (this%m_manager%x(j) == this%m_manager%x(j+1)) then
        rst = this%m_manager%y(j)
    else
        rst = this%m_manager%y(j) + ((x - this%m_manager%x(j)) / &
            (this%m_manager%x(j+1) - this%m_manager%x(j))) * &
            (this%m_manager%y(j+1) - this%m_manager(j))
    end if
end function

! ------------------------------------------------------------------------------
end module

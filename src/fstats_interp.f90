module fstats_interp
    use iso_fortran_env
    use fstats_errors
    use ferror
    implicit none
    private
    public :: interp_routine
    public :: base_interpolator
    public :: linear_interpolator

    type interpolation_manager
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
        procedure, public :: initialize => im_init
        procedure, public :: locate => im_locate
        procedure, public :: hunt => im_hunt
        procedure, public :: size => im_size
        procedure, public :: method => im_get_method
        procedure, public :: order => im_get_order
    end type

! ------------------------------------------------------------------------------
    type, abstract :: base_interpolator
        !! A base object for interpolation.
    contains
        procedure(interp_routine), public, deferred, pass :: interpolate_value
        procedure, public :: interpolate => bi_interp
    end type

    interface
        function interp_routine(this, x) result(rst)
            !! Performs a single interpolation.
            use iso_fortran_env, only : real64
            import base_interpolator
            class(base_interpolator), intent(inout) :: this
                !! The base_interpolator object.
            real(real64), intent(in) :: x
                !! The value at which to compute the interpolation.
            real(real64) :: rst
                !! The interpolated result.
        end function
    end interface

! ------------------------------------------------------------------------------
    type, extends(base_interpolator) :: linear_interpolator
        type(interpolation_manager), private :: m_manager
            !! An object to manage the interpolation process.
    contains
        procedure, public :: initialize => li_init
        procedure, public :: interpolate_value => li_raw_interp
    end type

contains
! ------------------------------------------------------------------------------
pure function is_monotonic(x) result(rst)
    !! Tests to see if an array is monotonic (either ascending or descending).
    real(real64), intent(in), dimension(:) :: x
        !! The array to test.
    logical :: rst
        !! Returns true if the array is monotonic; else, false.

    ! Local Variables
    integer(int32) :: i, n
    logical :: ascend

    ! Process
    n = size(x)
    ascend = x(n) > x(1)
    rst = .true.
    if (ascend) then
        do i = 2, n
            if (x(i) <= x(i-1)) then
                rst = .false.
                exit
            end if
        end do
    else
        do i = 2, n
            if (x(i) >= x(i-1)) then
                rst = .false.
                exit
            end if
        end do
    end if
end function

! ******************************************************************************
! INTERPOLATION_MANAGER
! ------------------------------------------------------------------------------
subroutine im_init(this, x, y, order, err)
    !! Initializes the interpolation object.
    class(interpolation_manager), intent(inout) :: this
        !! The interpolation_manager object.
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
        call report_array_size_error(errmgr, "im_init", "y", n, size(y))
        return
    end if

    if (order < 1) then
        call report_polynomial_order_error(errmgr, "im_init", order, 1)
        return
    end if

    if (.not.is_monotonic(x)) then
        call report_nonmonotonic_array_error(errmgr, "im_init", "x")
        return
    end if

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
    call report_memory_error(errmgr, "im_init", flag)
    return
end subroutine

! ------------------------------------------------------------------------------
pure function im_get_method(this) result(rst)
    !! Gets an integer representing the search method to utilize.
    class(interpolation_manager), intent(in) :: this
        !! The interpolation_manager object.
    integer(int32) :: rst
        !! The method identifier.
    rst = this%m_method
end function

! ------------------------------------------------------------------------------
pure function im_get_order(this) result(rst)
    !! Gets the polynomial order.
    class(interpolation_manager), intent(in) :: this
        !! The interpolation_manager object.
    integer(int32) :: rst
        !! The order.
    rst = this%m_order
end function

! ------------------------------------------------------------------------------
function im_locate(this, x) result(rst)
    !! Locates the index below the nearest entry to x.
    class(interpolation_manager), intent(inout) :: this
        !! The interpolation_manager object.
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
function im_hunt(this, x) result(rst)
    !! A search routine, similar to locate, but tailored towards searching
    !! when the initial location estimate is very near the desired point.  If
    !! not near, this method is rather inefficient; however, when near, this
    !! method is quite efficient when compared with locate.
    class(interpolation_manager), intent(inout) :: this
        !! The interpolation_manager object.
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
pure function im_size(this) result(rst)
    !! Gets the size of the stored raw data array.
    class(interpolation_manager), intent(in) :: this
        !! The interpolation_manager object.
    integer(int32) :: rst
        !! The size of the stored data.

    if (allocated(this%x)) then
        rst = size(this%x)
    else
        rst = 0
    end if
end function

! ******************************************************************************
! BASE_INTERPOLATOR
! ------------------------------------------------------------------------------
subroutine bi_interp(this, x, yi, err)
    !! Performs the interpolation.
    class(base_interpolator), intent(inout) :: this
        !! The base_interpolator object.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the x values at which to compute the
        !! interpolation.
    real(real64), intent(out), dimension(:) :: yi
        !! An N-element array  containing the interpolated data.
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.

    ! Local Variables
    integer(int32) :: i, n
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
        call report_array_size_error(errmgr, "bi_interp", "yi", n, size(yi))
        return
    end if

    ! Perform the interpolation
    do i = 1, n
        yi(i) = this%interpolate_value(x(i))
    end do
end subroutine

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
function li_raw_interp(this, x) result(rst)
    !! Raw, single-value interpolation routine.
    class(linear_interpolator), intent(inout) :: this
        !! The linear_interpolator object.
    real(real64), intent(in) :: x
        !! The value at which to compute the interpolation.
    real(real64) :: rst
        !! The interpolated value.

    ! Local Variables
    integer(int32) :: j

    ! Locate the proper index
    if (this%m_manager%method() == 1) then
        j = this%m_manager%hunt(x)
    else
        j = this%m_manager%locate(x)
    end if

    ! Perform the interpolation
    if (this%m_manager%x(j) == this%m_manager%x(j+1)) then
        rst = this%m_manager%y(j)
    else
        rst = this%m_manager%y(j) + ((x - this%m_manager%x(j)) / &
            (this%m_manager%x(j+1) - this%m_manager%x(j))) * &
            (this%m_manager%y(j+1) - this%m_manager%y(j))
    end if
end function

! ------------------------------------------------------------------------------
end module

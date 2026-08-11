module fstats_interp
    use iso_fortran_env
    use fstats_errors
    implicit none
    private
    public :: interp_routine
    public :: base_interpolator
    public :: linear_interpolator
    public :: polynomial_interpolator
    public :: spline_interpolator
    public :: SPLINE_QUADRATIC_OVER_INTERVAL
    public :: SPLINE_KNOWN_FIRST_DERIVATIVE
    public :: SPLINE_KNOWN_SECOND_DERIVATIVE
    public :: SPLINE_CONTINUOUS_THIRD_DERIVATIVE
    public :: hermite_interpolator

        
    integer(int32), parameter :: SPLINE_QUADRATIC_OVER_INTERVAL = 1000
        !! Indicates that the spline is quadratic over the interval under
        !! consideration (beginning or ending interval).
    integer(int32), parameter :: SPLINE_KNOWN_FIRST_DERIVATIVE = 1001
        !! Indicates a known first derivative at either the beginning or ending
        !! point.
    integer(int32), parameter :: SPLINE_KNOWN_SECOND_DERIVATIVE = 1002
        !! Indicates a known second derivative at either the beginning or ending
        !! point.
    integer(int32), parameter :: SPLINE_CONTINUOUS_THIRD_DERIVATIVE = 1003
        !! Indicates a continuous third derivative at either the beginning or 
        !! ending point.

! ------------------------------------------------------------------------------
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
        !! Defines a type meant for performing piecewise linear interpolation.
        type(interpolation_manager), private :: m_manager
            !! An object to manage the interpolation process.
    contains
        procedure, public :: initialize => li_init
        procedure, public :: interpolate_value => li_raw_interp
    end type

! ------------------------------------------------------------------------------
    type, extends(base_interpolator) :: polynomial_interpolator
        !! Defines a type meant for performing piecewise polynomial
        !! interpolation.
        type(interpolation_manager), private :: m_manager
            !! An object to manage the interpolation process.
        real(real64), private, allocatable, dimension(:) :: m_c
            !! Workspace array.
        real(real64), private, allocatable, dimension(:) :: m_d
            !! Workspace array.
    contains
        procedure, public :: initialize => pi_init
        procedure, public :: interpolate_value => pi_raw_interp
    end type

! ------------------------------------------------------------------------------
    type, extends(base_interpolator) :: spline_interpolator
        !! Defines a type meant for performing cubic spline interpolation.
        type(interpolation_manager), private :: m_manager
            !! An object to manage the interpolation process.
        real(real64), private, allocatable, dimension(:) :: m_ypp
            !! A workspace array.
    contains
        procedure, public :: initialize => si_init
        procedure, public :: interpolate_value => si_raw_interp
        procedure, private :: compute_second_derivatives => si_second_diff
    end type

! ------------------------------------------------------------------------------
    type, extends(base_interpolator) :: hermite_interpolator
        !! Defines a type meant for performing Hermite-type interpolation.
        !! The interpolating polynomial constructed by this object is a global
        !! polynomial, not a piecewise polynomial.  Given N data points, the
        !! polynoial will be of degree 2 * N - 1.  As N increases, the
        !! interpolating polynomial may be liable to oscillations that do
        !! not properly represent the data.  For large data sets, a piecewise
        !! polynomial approach is recommended.  See either the
        !! polynomial_interpolator or spline_interpolator types.
        !!
        !! This implementation is a modification of the HERMITE library
        !! which can be found 
        !! [here](https://people.math.sc.edu/Burkardt/f_src/hermite/hermite.html).
        real(real64), private, allocatable, dimension(:) :: m_x
            !! The x-coordinate raw data.
        real(real64), private, allocatable, dimension(:) :: m_y
            !! The y-coordinate raw data.
        real(real64), private, allocatable, dimension(:) :: m_yp
            !! An N-element array containing the first derivatives.
        real(real64), private, allocatable, dimension(:) :: m_xd
            !! A 2*N-element array containing the abscissas for the divided 
            !! difference table.
        real(real64), private, allocatable, dimension(:) :: m_yd
            !! A 2*N-element array containing the divided difference table.
        real(real64), private, allocatable, dimension(:) :: m_xdp
            !! A 2*N-1-element array containing the derivatives for the
            !! abscissas for the divided difference table.
        real(real64), private, allocatable, dimension(:) :: m_ydp
            !! A 2*N-1-element array containing the derivatives for the
            !! divided difference table.
        real(real64), private, allocatable, dimension(:) :: m_xwork
            !! A workspace array.
        real(real64), private, allocatable, dimension(:) :: m_ywork
            !! A workspace array.
    contains
        procedure, private :: compute_ddf_derivatives => hi_dif_deriv
        procedure, private :: set_up_hermite_interpolant => hi_set_up_table
        procedure, public :: initialize => hi_init
        procedure, public :: interpolate_value => hi_raw_interp
        procedure, public :: interpolate => hi_interp
        procedure, public :: interpolate_with_derivative => hi_interp_all
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
subroutine im_init(this, x, y, order)
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

    ! Local Variables
    integer(int32) :: i, n
    
    ! Initialization
    n = size(x)

    ! Input Checking
    if (size(y) /= n) error stop FS_ARRAY_SIZE_ERROR

    if (order < 1) error stop 4

    if (.not.is_monotonic(x)) error stop FS_NONMONOTONIC_ARRAY_ERROR

    ! Process
    this%m_order = order
    this%m_nearestIndex = 1
    this%m_indexDelta = 1
    this%m_method = 0

    if (allocated(this%x)) deallocate(this%x)
    if (allocated(this%y)) deallocate(this%y)
    allocate(this%x(n), this%y(n))
    do i = 1, n
        this%x(i) = x(i)
        this%y(i) = y(i)
    end do
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
subroutine bi_interp(this, x, yi)
    !! Performs the interpolation.
    class(base_interpolator), intent(inout) :: this
        !! The base_interpolator object.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the x values at which to compute the
        !! interpolation.
    real(real64), intent(out), dimension(:) :: yi
        !! An N-element array  containing the interpolated data.

    ! Local Variables
    integer(int32) :: i, n
    
    ! Initialization
    n = size(x)

    ! Input Checking
    if (size(yi) /= n) error stop FS_ARRAY_SIZE_ERROR

    ! Perform the interpolation
    do i = 1, n
        yi(i) = this%interpolate_value(x(i))
    end do
end subroutine

! ******************************************************************************
! LINEAR_INTERPOLATOR
! ------------------------------------------------------------------------------
subroutine li_init(this, x, y)
    !! Initializes the interpolation object.
    class(linear_interpolator), intent(inout) :: this
        !! The linear_interpolator object.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the x-coordinate data in either 
        !! monotonically increasing or decreasing order.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the y-coordinate data.

    ! Initialize the managing object
    call this%m_manager%initialize(x, y, 1)
end subroutine

! ------------------------------------------------------------------------------
function li_raw_interp(this, x) result(rst)
    !! Interpolates a single value.
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

! ******************************************************************************
! POLYNOMIAL_INTERPOLATOR
! ------------------------------------------------------------------------------
function pi_raw_interp(this, x) result(rst)
    !! Interpolates a single value.
    class(polynomial_interpolator), intent(inout) :: this
        !! The polynomial_interpolator object.
    real(real64), intent(in) :: x
        !! The value at which to compute the interpolation.
    real(real64) :: rst
        !! The interpolated value.

    ! Local Variables
    integer(int32) :: i, ind, m, ns, mm, jl, jlo
    real(real64) :: den, dif, dift, ho, hp, w, dy

    ! Initialization
    if (this%m_manager%method() == 1) then
        jlo = this%m_manager%hunt(x)
    else
        jlo = this%m_manager%locate(x)
    end if
    mm = this%m_manager%order() + 1
    ns = 1
    if (jlo == 1) then
        jl = 1
    else
        jl = jlo - 1
    end if
    dif = abs(x - this%m_manager%x(jl))

    ! Process
    do i = 1, mm
        ind = jl + i - 1
        dift = abs(x - this%m_manager%x(ind))
        if (dift < dif) then
            ns = i
            dif = dift
        end if
        this%m_c(i) = this%m_manager%y(ind)
        this%m_d(i) = this%m_manager%y(ind)
    end do

    rst = this%m_manager%y(jl + ns - 1)
    ns = ns - 1

    do m = 1, mm - 1
        do i = 1, mm - m
            ind = jl + i - 1
            ho = this%m_manager%x(ind) - x
            hp = this%m_manager%x(ind + m) - x
            w = this%m_c(i + 1) - this%m_d(i)
            den = ho - hp
            den = w / den
            this%m_d(i) = hp * den
            this%m_c(i) = ho * den
        end do
        if (2 * ns < mm - m) then
            dy = this%m_c(ns + 1)
        else
            dy = this%m_d(ns)
            ns = ns - 1
        end if
        rst = rst + dy
    end do
end function

! ------------------------------------------------------------------------------
subroutine pi_init(this, order, x, y)
    !! Initializes the interpolation object.
    class(polynomial_interpolator), intent(inout) :: this
        !! The polynomial_interpolator object.
    integer(int32), intent(in) :: order
        !! The polynomial order.  This value must be at least 1.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the x-coordinate data in either
        !! monotonically increasing or decreasing order.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the y-coordinate data.

    ! Local Variables
    integer(int32) :: m
    
    ! Initialization
    call this%m_manager%initialize(x, y, order)

    ! Allocate necessary workspace memory
    if (allocated(this%m_c)) deallocate(this%m_c)
    if (allocated(this%m_d)) deallocate(this%m_d)
    m = order + 1
    allocate(this%m_c(m), this%m_d(m))
end subroutine

! ******************************************************************************
! SPLINE_INTERPOLATOR
! ------------------------------------------------------------------------------
function si_raw_interp(this, x) result(rst)
    !! Interpolates a single value.
    class(spline_interpolator), intent(inout) :: this
        !! The spline_interpolator object.
    real(real64), intent(in) :: x
        !! The value at which to compute the interpolation.
    real(real64) :: rst
        !! The interpolated value.

    ! Local Variables
    integer(int32) :: jlo, right
    real(real64) :: dt, h, yr, yl

    ! Initialization
    if (this%m_manager%method() == 1) then
        jlo = this%m_manager%hunt(x)
    else
        jlo = this%m_manager%locate(x)
    end if
    right = jlo + 1
    dt = x - this%m_manager%x(jlo)
    h = this%m_manager%x(right) - this%m_manager%x(jlo)

    ! Process
    yr = this%m_manager%y(right)
    yl = this%m_manager%y(jlo)
    rst = yl + dt * ((yr - yl) / h - (this%m_ypp(right) / 6.0d0 + &
        this%m_ypp(jlo) / 3.0d0) * h + dt * (0.5d0 * this%m_ypp(jlo) + &
        dt * ((this%m_ypp(right) - this%m_ypp(jlo)) / (6.0d0 * h))))
end function

! ------------------------------------------------------------------------------
! [Spline Library](http://people.sc.fsu.edu/~jburkardt/f77_src/spline/spline.html)
subroutine penta_solve(a1, a2, a3, a4, a5, b, x)
    !! Solves a pentadiagonal system of linear equations.  A
    !! pentadiagonal matrix is all zeros with the exception of the diagonal,
    !! and the two immediate sub and super-diagonals.  The entries of row I
    !! are stored as follows:
    !!      A(I,I-2) -> A1(I)
    !!      A(I,I-1) -> A2(I)
    !!      A(I,I) -> A3(I)
    !!      A(I,I+1) -> A4(I)
    !!      A(I,I+2) -> A5(I)
    real(real64), intent(in), dimension(:) :: a1
        !! An N-element array as defined above.
    real(real64), intent(inout), dimension(:) :: a2
        !! An N-element array as defined above.  This array is
        !! overwritten by this routine during the solution process.
    real(real64), intent(inout), dimension(:) :: a3
        !! An N-element array as defined above.  This array is
        !! overwritten by this routine during the solution process.
    real(real64), intent(inout), dimension(:) :: a4
        !! An N-element array as defined above.  This array is
        !! overwritten by this routine during the solution process.
    real(real64), intent(in), dimension(:) :: a5
        !! An N-element array as defined above.
    real(real64), intent(inout), dimension(:) :: b
        !! An N-element array containing the right-hand-side.  This array is 
        !! overwritten by this routine during the solution process.
    real(real64), intent(out), dimension(:) :: x
        !! An N-element array that, on output, contains the solution to the 
        !! linear system.

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: xmult

    ! Initialization
    n = size(a1)

    ! Process
    do i = 2, n - 1
        xmult = a2(i) / a3(i - 1)
        a3(i) = a3(i) - xmult * a4(i - 1)
        a4(i) = a4(i) - xmult * a5(i - 1)
        b(i) = b(i) - xmult * b(i - 1)
        xmult = a1(i + 1) - xmult * a4(i - 1)
        a2(i + 1) = a2(i + 1) - xmult * a4(i - 1)
        a3(i + 1) = a3(i + 1) - xmult * a5(i - 1)
        b(i + 1) = b(i + 1) - xmult * b(i - 1)
    end do

    xmult = a2(n) / a3(n - 1)
    a3(n) = a3(n) - xmult * a4(n - 1)
    x(n) = (b(n) - xmult * b(n - 1)) / a3(n)
    x(n - 1) = (b(n - 1) - a4(n - 1) * x(n)) / a3(n - 1)
    do i = n - 2, 1, -1
        x(i) = (b(i) - a4(i) * x(i + 1) - a5(i) * x(i + 2)) / a3(i)
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine si_init(this, x, y, ibcbeg, ybcbeg, ibcend, ybcend)
    !! Initializes the interpolation object.
    class(spline_interpolator), intent(inout) :: this
        !! The spline_interpolator object.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the x-coordinate data in either
        !! monotonically increasing or decreasing order.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the y-coordinate data.
    integer(int32), intent(in), optional :: ibcbeg
        !! An optional input that defines the nature of the boundary condition 
        !! at the beginning of the spline.  If no parameter, or an invalid 
        !! parameter, is specified, the default condition
        !!  (SPLINE_QUADRATIC_OVER_INTERVAL) is used.
        !!
        !!  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
        !!      initial interval.  No value is required for ybcbeg.
        !!
        !!  - SPLINE_KNOWN_FIRST_DERIVATIVE: The spline's first derivative at 
        !!      its initial point is provided in ybcbeg.
        !!
        !!  - SPLINE_KNOWN_SECOND_DERIVATIVE: The spline's second derivative at 
        !!      its initial point is provided in ybcbeg.
        !!
        !!  - SPLINE_CONTINUOUS_THIRD_DERIVATIVE: The third derivative is 
        !!      continuous at x(2).  No value is required for ybcbeg.
    real(real64), intent(in), optional :: ybcbeg
        !! If needed, the value of the initial point boundary condition.  If 
        !! needed, but not supplied, a default value of zero will be used.
    integer(int32), intent(in), optional :: ibcend
        !! An optional input that defines the nature of the boundary condition 
        !! at the end of the spline.  If no parameter, or an invalid 
        !! parameter, is specified, the default condition
        !!  (SPLINE_QUADRATIC_OVER_INTERVAL) is used.
        !!
        !!  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
        !!      final interval.  No value is required for ybcend.
        !!
        !!  - SPLINE_KNOWN_FIRST_DERIVATIVE: The spline's first derivative at 
        !!      its final point is provided in ybcend.
        !!
        !!  - SPLINE_KNOWN_SECOND_DERIVATIVE: The spline's second derivative at 
        !!      its final point is provided in ybcend.
        !!
        !!  - SPLINE_CONTINUOUS_THIRD_DERIVATIVE: The third derivative is 
        !!      continuous at x(n-1).  No value is required for ybcend.
    real(real64), intent(in), optional :: ybcend
        !! If needed, the value of the final point boundary condition.  If 
        !! needed, but not supplied, a default value of zero will be used.

    ! Local Variables
    integer(int32) :: ibeg, iend, n
    real(real64) :: ybeg, yend
    
    ! Initialization
    n = size(x)
    ibeg = SPLINE_QUADRATIC_OVER_INTERVAL
    iend = SPLINE_QUADRATIC_OVER_INTERVAL
    ybeg = 0.0d0
    yend = 0.0d0
    if (present(ibcbeg)) ibeg = ibcbeg
    if (present(ibcend)) iend = ibcend
    if (present(ybcbeg)) ybeg = ybcbeg
    if (present(ybcend)) yend = ybcend

    ! Input Checking
    if (size(y) /= n) error stop FS_ARRAY_SIZE_ERROR
    if (ibeg /= SPLINE_CONTINUOUS_THIRD_DERIVATIVE .and. &
        ibeg /= SPLINE_KNOWN_SECOND_DERIVATIVE .and. &
        ibeg /= SPLINE_KNOWN_FIRST_DERIVATIVE .and. &
        ibeg /= SPLINE_QUADRATIC_OVER_INTERVAL) &
    then
        ibeg = SPLINE_QUADRATIC_OVER_INTERVAL
    end if
    if (iend /= SPLINE_CONTINUOUS_THIRD_DERIVATIVE .and. &
        iend /= SPLINE_KNOWN_SECOND_DERIVATIVE .and. &
        iend /= SPLINE_KNOWN_FIRST_DERIVATIVE .and. &
        iend /= SPLINE_QUADRATIC_OVER_INTERVAL) &
    then
        iend = SPLINE_QUADRATIC_OVER_INTERVAL
    end if
    call this%m_manager%initialize(x, y, 3)

    ! Evaluate the 2nd derivatives
    call this%compute_second_derivatives(ibeg, ybeg, iend, yend)
end subroutine

! ------------------------------------------------------------------------------
! http://people.sc.fsu.edu/~jburkardt/f77_src/spline/spline.html
! ROUTINE: SPLINE_CUBIC_SET
subroutine si_second_diff(this, ibcbeg, ybcbeg, ibcend, ybcend)
    !! Computes the second derivative terms for the cubic-spline model.
    class(spline_interpolator), intent(inout) :: this
        !! The spline_interpolator object.
    integer(int32), intent(in) :: ibcbeg
        !! An input that defines the nature of the boundary condition 
        !! at the beginning of the spline.
        !!
        !!  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
        !!      initial interval.
        !!
        !!  - SPLINE_KNOWN_FIRST_DERIVATIVE: The spline's first derivative at 
        !!      its initial point is provided in ybcbeg.
        !!
        !!  - SPLINE_KNOWN_SECOND_DERIVATIVE: The spline's second derivative at 
        !!      its initial point is provided in ybcbeg.
        !!
        !!  - SPLINE_CONTINUOUS_THIRD_DERIVATIVE: The third derivative is 
        !!      continuous at x(2).
    real(real64), intent(in) :: ybcbeg
        !! The value of the initial point boundary condition.
    integer(int32), intent(in) :: ibcend
        !! An input that defines the nature of the boundary condition 
        !! at the end of the spline.
        !!
        !!  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
        !!      final interval.
        !!
        !!  - SPLINE_KNOWN_FIRST_DERIVATIVE: The spline's first derivative at 
        !!      its final point is provided in ybcbeg.
        !!
        !!  - SPLINE_KNOWN_SECOND_DERIVATIVE: The spline's second derivative at 
        !!      its final point is provided in ybcbeg.
        !!
        !!  - SPLINE_CONTINUOUS_THIRD_DERIVATIVE: The third derivative is 
        !!      continuous at x(n-1).
    real(real64), intent(in) :: ybcend
        !! The value of the final point boundary condition.

    ! Local Variables
    integer(int32) :: i, n
    real(real64), allocatable ,dimension(:) :: a1, a2, a3, a4, a5, b
    
    ! Initialization
    n = this%m_manager%size()

    ! Local Memory Allocations
    if (allocated(this%m_ypp)) deallocate(this%m_ypp)
    allocate(this%m_ypp(n), a1(n), a2(n), a3(n), a4(n), a5(n), b(n), source = 0.0d0)

    ! Set up the first equation
    select case (ibcbeg)
    case (SPLINE_QUADRATIC_OVER_INTERVAL)
        b(1) = 0.0d0
        a3(1) = 1.0d0
        a4(1) = -1.0d0
    case (SPLINE_KNOWN_FIRST_DERIVATIVE)
        b(1) = (this%m_manager%y(2) - this%m_manager%y(1)) / &
            (this%m_manager%x(2) - this%m_manager%x(1)) - ybcbeg
        a3(1) = (this%m_manager%x(2) - this%m_manager%x(1)) / 3.0d0
        a4(1) = (this%m_manager%x(2) - this%m_manager%x(1)) / 6.0d0
    case (SPLINE_KNOWN_SECOND_DERIVATIVE)
        b(1) = ybcbeg
        a3(1) = 1.0d0
        a4(1) = 0.0d0
    case default
        b(1) = 0.0d0
        a3(1) = this%m_manager%x(2) - this%m_manager%x(3)
        a4(1) = this%m_manager%x(3) - this%m_manager%x(1)
        a5(1) = this%m_manager%x(1) - this%m_manager%x(2)
    end select

    ! Set up the intermediate equations
    do i = 2, n - 1
        b(i) = (this%m_manager%y(i+1) - this%m_manager%y(i)) / &
            (this%m_manager%x(i+1) - this%m_manager%x(i)) - &
            (this%m_manager%y(i) - this%m_manager%y(i-1)) / &
            (this%m_manager%x(i) - this%m_manager%x(i-1))
        a2(i) = (this%m_manager%x(i+1) - this%m_manager%x(i)) / 6.0d0
        a3(i) = (this%m_manager%x(i+1) - this%m_manager%x(i-1)) / 3.0d0
        a4(i) = (this%m_manager%x(i) - this%m_manager%x(i-1)) / 6.0d0
    end do

    ! Set up the last equation
    select case (ibcend)
    case (SPLINE_QUADRATIC_OVER_INTERVAL)
        b(n) = 0.0d0
        a2(n) = -1.0d0
        a3(n) = 1.0d0
    case (SPLINE_KNOWN_FIRST_DERIVATIVE)
        b(n) = ybcend - (this%m_manager%y(n) - this%m_manager%y(n-1)) / &
            (this%m_manager%x(n) - this%m_manager%x(n-1))
        a2(n) = (this%m_manager%x(n) - this%m_manager%x(n-1)) / 6.0d0
        a3(n) = (this%m_manager%x(n) - this%m_manager%x(n-1)) / 3.0d0
    case (SPLINE_KNOWN_SECOND_DERIVATIVE)
        b(n) = ybcend
        a2(n) = 0.0d0
        a3(n) = 1.0d0
    case (SPLINE_CONTINUOUS_THIRD_DERIVATIVE)
        b(n) = 0.0d0
        a1(n) = this%m_manager%x(n-1) - this%m_manager%x(n)
        a2(n) = this%m_manager%x(n) - this%m_manager%x(n-2)
        a3(n) = this%m_manager%x(n-2) - this%m_manager%x(n-1)
    case default
        b(n) = 0.0d0
        a2(n) = -1.0d0
        a3(n) = 1.0d0
    end select

    ! Define the 2nd derivative
    if (n == 2 .and. ibcbeg == SPLINE_QUADRATIC_OVER_INTERVAL .and. &
        ibcend == SPLINE_QUADRATIC_OVER_INTERVAL) &
    then
        this%m_ypp(1) = 0.0d0
        this%m_ypp(2) = 0.0d0
    else
        call penta_solve(a1, a2, a3, a4, a5, b, this%m_ypp)
    end if
end subroutine

! ******************************************************************************
! HERMITE_INTERPOLATOR
! ------------------------------------------------------------------------------
! REF: https://people.math.sc.edu/Burkardt/f_src/hermite/hermite.html
subroutine hi_dif_deriv(this)
    !! Computes the derivatives of a polynomial in divided difference form.
    class(hermite_interpolator), intent(inout) :: this
        !! The hermite_interpolator object.

    ! Local Variables
    integer(int32) :: i, nd, ndp
    
    ! Initialization
    nd = size(this%m_xd)
    ndp = nd - 1

    ! Memory Allocations
    if (allocated(this%m_xwork)) deallocate(this%m_xwork)
    if (allocated(this%m_ywork)) deallocate(this%m_ywork)
    if (allocated(this%m_xdp)) deallocate(this%m_xdp)
    if (allocated(this%m_ydp)) deallocate(this%m_ydp)
    allocate(this%m_xwork(nd), this%m_ywork(nd), this%m_xdp(ndp), &
        this%m_ydp(ndp), source = 0.0d0)

    ! Create a copy of the difference table and shift the abscissas to zero
    this%m_xwork = this%m_xd
    this%m_ywork = this%m_yd
    call dif_shift_zero(this%m_xwork, this%m_ywork)

    ! Construct the derivative
    this%m_ydp = (/ (i * this%m_ywork(i), i = 1, ndp) /)
end subroutine

! ------------------------------------------------------------------------------
subroutine hi_set_up_table(this)
    !! Sets up the divided difference table from the raw data.
    class(hermite_interpolator), intent(inout) :: this
        !! The hermite_interpolator object.

    ! Local Variables
    integer(int32) :: i, j, n, nd, ndp
    
    ! Initialization
    n = size(this%m_x)
    nd = 2 * n
    ndp = nd - 1

    ! Memory Allocations
    allocate(this%m_xd(nd), this%m_yd(nd))

    ! Copy data
    this%m_xd(1:nd-1:2) = this%m_x(1:n)
    this%m_xd(2:nd:2) = this%m_x(1:n)

    ! Start differencing
    this%m_yd(1) = this%m_y(1)
    this%m_yd(3:nd-1:2) = (this%m_y(2:n) - this%m_y(1:n-1)) / &
        (this%m_x(2:n) - this%m_x(1:n-1))
    this%m_yd(2:nd:2) = this%m_yp(1:n)

    ! Carry out the remaining steps in the usual way
    do i = 3, nd
        do j = nd, i, -1
            this%m_yd(j) = (this%m_yd(j) - this%m_yd(j-1)) / &
                (this%m_xd(j) - this%m_xd(j+1-i))
        end do
    end do

    ! Compute the difference table for the derivative
    call this%compute_ddf_derivatives()
end subroutine

! ------------------------------------------------------------------------------
subroutine hi_init(this, x, y, yp)
    !! Initializes the interpolation object.
    class(hermite_interpolator), intent(inout) :: this
        !! The hermite_interpolator object.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the x-coordinate data in either
        !! monotonically increasing or decreasing order.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the y-coordinate data.
    real(real64), intent(in), dimension(:) :: yp
        !! An N-element array containing the first derivative of the data.

    ! Local Variables
    integer(int32) :: n
    
    ! Initialization
    n = size(x)

    ! Input Checking
    if (size(y) /= n) error stop 3
    if (size(yp) /= n) error stop 4

    ! Allocate memory
    if (allocated(this%m_x)) deallocate(this%m_x)
    if (allocated(this%m_y)) deallocate(this%m_y)
    if (allocated(this%m_yp)) deallocate(this%m_yp)
    allocate(this%m_x(n), source = x)
    allocate(this%m_y(n), source = y)
    allocate(this%m_yp(n), source = yp)

    ! Set up the table
    call this%set_up_hermite_interpolant()
end subroutine

! ------------------------------------------------------------------------------
function hi_raw_interp(this, x) result(rst)
    !! Interpolates a single value.
    class(hermite_interpolator), intent(inout) :: this
        !! The hermite_interpolator object.
    real(real64), intent(in) :: x
        !! The value at which to compute the interpolation.
    real(real64) :: rst
        !! The interpolated value.

    ! Local Variables
    real(real64) :: yi(1)

    ! Process
    call hi_interp(this, [x], yi)
    rst = yi(1)
end function

! ------------------------------------------------------------------------------
subroutine hi_interp_all(this, x, yi, ypi)
    !! Performs the interpolation.
    class(hermite_interpolator), intent(inout) :: this
        !! The hermite_interpolator object.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the x values at which to compute the
        !! interpolation.
    real(real64), intent(out), dimension(:) :: yi
        !! An N-element array  containing the interpolated data.
    real(real64), intent(out), optional, dimension(:) :: ypi
        !! An N-element array containing the interpolated first derivative
        !! data, if supplied.

    ! Local Variables
    integer(int32) :: i, nd, nv, ndp
    
    ! Initialization
    nd = size(this%m_xd)
    ndp = nd - 1
    nv = size(x)

    ! Ensure everything is initialized
    if (.not.allocated(this%M_x)) error stop FS_UNINITIALIZED_OBJECT_ERROR

    ! Input Check
    if (size(yi) /= nv) error stop 3
    if (present(ypi)) then
        if (size(ypi) /= nv) error stop 4
    end if

    ! Process
    yi(1:nv) = this%m_yd(nd)
    do i = nd - 1, 1, -1
        yi(1:nv) = this%m_yd(i) + (x(1:nv) - this%m_xd(i)) * yi(1:nv)
    end do

    ! Derivative?
    if (present(ypi)) then
        ypi(1:nv) = this%m_ydp(ndp)
        do i = ndp - 1, 1, -1
            ypi(1:nv) = this%m_ydp(i) + (x(1:nv) - this%m_xdp(i)) * ypi(1:nv)
        end do
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine hi_interp(this, x, yi)
    !! Performs the interpolation.
    class(hermite_interpolator), intent(inout) :: this
        !! The hermite_interpolator object.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the x values at which to compute the
        !! interpolation.
    real(real64), intent(out), dimension(:) :: yi
        !! An N-element array  containing the interpolated data.

    ! Process
    call hi_interp_all(this, x, yi)
end subroutine

! ------------------------------------------------------------------------------
subroutine dif_shift_x(xd, yd, xv)
    !! Replaces one abcissa of a divided difference table.
    real(real64), intent(inout), dimension(:) :: xd
        !! The x-values used in the representation of the divided difference
        !! polynomial.  On output, the last entry of this array has been dropped
        !! and the other entries have shifted up one index.  The value xv has 
        !! been inserted at the beginning of the array.
    real(real64), intent(inout), dimension(:) :: yd
        !! The divided difference coefficients corresponding to xd.  On output,
        !! this array has beed adjusted accordingly.
    real(real64), intent(in) :: xv
        !! A new x-value which is to be used in the representation of the
        !! polynomial.

    ! Local Variables
    integer(int32) :: i, nd

    ! Recompute the divided difference coefficients
    nd = size(xd)
    do i = nd - 1, 1, -1
        yd(i) = yd(i) + (xv - xd(i)) * yd(i + 1)
    end do

    ! Shift XD
    xd(2:nd) = xd(1:nd-1)
    xd(1) = xv
end subroutine

! ------------------------------------------------------------------------------
subroutine dif_shift_zero(xd, yd)
    !! Shifts a divided difference table so all abscissas are zero.
    real(real64), intent(inout), dimension(:) :: xd
        !! The x-values that correspond to the divided difference table.  On
        !! output, this array contains only zeros.
    real(real64), intent(inout), dimension(:) :: yd
        !! The divided difference table.  On putput, this array is also the 
        !! coefficient array for the standard representation of the polynomial.

    ! Local Variables
    integer(int32) :: i, nd
    real(real64) :: xv

    ! Process
    nd = size(xd)
    xv = 0.0d0
    do i = 1, nd
        call dif_shift_x(xd, yd, xv)
    end do
end subroutine

! ------------------------------------------------------------------------------
end module

module fstats_mcmc
    use iso_fortran_env
    use fstats_distributions
    use fstats_sampling
    use ferror
    use fstats_errors
    use fstats_descriptive_statistics
    use fstats_types
    use fstats_regression
    implicit none
    private
    public :: metropolis_hastings

    type metropolis_hastings
        !! An implementation of the Metropolis-Hastings algorithm for the
        !! generation of a Markov chain.
        integer(int32), private :: initial_iteration_estimate = 10000
            !! An initial estimate at the number of allowed iterations.
        integer(int32), private :: m_bufferSize = 0
            !! The actual number of states (# of used rows) in the buffer.
        real(real64), private, allocatable, dimension(:,:) :: m_buffer
            !! The buffer where each new state is stored as a row in the matrix.
        integer(int32), private :: m_numVars = 0
            !! The number of state variables.
    contains
        procedure, public :: get_state_variable_count => mh_get_nvars
        procedure, public :: get_chain_length => mh_get_chain_length
        procedure, public :: push_new_state => mh_push
        procedure, public :: get_chain => mh_get_chain
        procedure, private :: resize_buffer => mh_resize_buffer
        procedure, private :: get_buffer_length => mh_get_buffer_length
    end type

contains
! ------------------------------------------------------------------------------
pure function mh_get_nvars(this) result(rst)
    !! Gets the number of state variables.
    class(metropolis_hastings), intent(in) :: this
        !! The metropolis_hastings object.
    integer(int32) :: rst
        !! The number of state variables.

    rst = this%m_numVars
end function

! ------------------------------------------------------------------------------
pure function mh_get_chain_length(this) result(rst)
    !! Gets the length of the chain (number of stored state variables).
    class(metropolis_hastings), intent(in) :: this
        !! The metropolis_hastings object.
    integer(int32) :: rst
        !! The chain length.

    rst = this%m_bufferSize
end function

! ------------------------------------------------------------------------------
subroutine mh_resize_buffer(this, err)
    !! Resizes the buffer to accept more states.
    class(metropolis_hastings), intent(inout) :: this
        !! The metropolis_hastings object.
    class(errors), intent(inout), optional, target :: err
        !! The error handling object.

    ! Local Variables
    integer(int32) :: m, n, flag, mOld
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    real(real64), allocatable, dimension(:,:) :: copy
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    m = this%initial_iteration_estimate
    n = this%get_state_variable_count()

    ! Is this the first time?
    if (.not.allocated(this%m_buffer)) then
        allocate(this%m_buffer(m, n), stat = flag)
        if (flag /= 0) go to 10
        return
    end if

    ! If we're here, then we need to create a copy and go from there
    m = size(this%m_buffer, 1)
    mOld = m
    allocate(copy(m, n), stat = flag, source = this%m_buffer)
    if (flag /= 0) go to 10
    deallocate(this%m_buffer)
    m = m + this%initial_iteration_estimate
    allocate(this%m_buffer(m, n), stat = flag)
    if (flag /= 0) go to 10
    this%m_buffer(1:mOld,:) = copy
    deallocate(copy)

    ! End
    return

    ! Memory Error Handling
10  continue
    call report_memory_error(errmgr, "mh_resize_buffer", flag)
    return
end subroutine

! ------------------------------------------------------------------------------
pure function mh_get_buffer_length(this) result(rst)
    !! Gets the actual length of the buffer.  This value will likely exceed the
    !! actual number of items in the chain.
    class(metropolis_hastings), intent(in) :: this
        !! The metropolis_hastings object.
    integer(int32) :: rst
        !! The actual buffer length.

    if (allocated(this%m_buffer)) then
        rst = size(this%m_buffer, 1)
    else
        rst = 0
    end if
end function

! ------------------------------------------------------------------------------
subroutine mh_push(this, x, err)
    !! Pushes a new set of state variables onto the buffer.
    class(metropolis_hastings), intent(inout) :: this
        !! The metropolis_hastings object.
    real(real64), intent(in), dimension(:) :: x
        !! The new N-element state array.
    class(errors), intent(inout), optional, target :: err
        !! The error handling object.

    ! Local Variables
    integer(int32) :: n, n1, nbuffer, nvars
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = this%get_chain_length()
    n1 = n + 1
    nbuffer = this%get_buffer_length()
    nvars = size(x)

    ! If this is the first time, ensure the collection is initialized
    if (this%get_state_variable_count() == 0) then
        this%m_numVars = nvars
    end if

    ! Input Checking
    if (nvars /= this%get_state_variable_count()) then
        call report_array_size_error(errmgr, "mh_push", "x", &
            this%get_state_variable_count(), nvars)
        return
    end if

    ! Resize the buffer, if necessary
    if (n == 0 .or. n == nbuffer) then
        call this%resize_buffer(errmgr)
        if (errmgr%has_error_occurred()) return
    end if

    ! Store the new state
    this%m_buffer(n1,:) = x
    this%m_bufferSize = this%m_bufferSize + 1
end subroutine

! ------------------------------------------------------------------------------
function mh_get_chain(this, err) result(rst)
    !! Gets a copy of the stored Markov chain.
    class(metropolis_hastings), intent(in) :: this
        !! The metropolis_hastings object.
    class(errors), intent(inout), optional, target :: err
        !! The error handling object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting chain with each parameter represented by a column.

    ! Local Variables
    integer(int32) :: npts, nvar, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    npts = this%get_chain_length()
    nvar = this%get_state_variable_count()

    ! Process
    allocate(rst(npts, nvar), stat = flag, &
        source = this%m_buffer(1:npts,1:nvar))
    if (flag /= 0) then
        call report_memory_error(errmgr, "mh_get_chain", flag)
        return
    end if
end function

! ------------------------------------------------------------------------------
end module
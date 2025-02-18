module fstats_mcmc
    use iso_fortran_env
    use fstats_distributions
    use fstats_sampling
    use ferror
    use fstats_errors
    use fstats_descriptive_statistics
    use fstats_types
    use fstats_regression
    use fstats_sampling
    implicit none
    private
    public :: metropolis_hastings

    type metropolis_hastings
        !! An implementation of the Metropolis-Hastings algorithm for the
        !! generation of a Markov chain.  This is a default implementation
        !! that allows sampling of normally distributed posterior distributions
        !! centered on zero with unit standard deviations.  Proposals are
        !! generated from a multivariate normal distribution with an identity
        !! covariance matrix and centered on zero.  To alter these sampling
        !! and target distributions simply create a new class inheriting from 
        !! this class and override the appropriate routines.
        integer(int32), private :: initial_iteration_estimate = 10000
            !! An initial estimate at the number of allowed iterations.
        integer(int32), private :: m_bufferSize = 0
            !! The actual number of states (# of used rows) in the buffer.
        real(real64), private, allocatable, dimension(:,:) :: m_buffer
            !! The buffer where each new state is stored as a row in the matrix.
        integer(int32), private :: m_numVars = 0
            !! The number of state variables.
        type(multivariate_normal_distribution) :: m_propDist
            !! The proposal distribution from which to draw samples.
    contains
        procedure, public :: get_state_variable_count => mh_get_nvars
        procedure, public :: get_chain_length => mh_get_chain_length
        procedure, public :: push_new_state => mh_push
        procedure, public :: get_chain => mh_get_chain
        procedure, public :: generate_proposal => mh_proposal
        procedure, public :: compute_hastings_ratio => mh_hastings_ratio
        procedure, public :: target_distribution => mh_target
        procedure, public :: sample => mh_sample
        procedure, public :: reset => mh_clear_chain
        procedure, public :: on_acceptance => mh_on_success
        procedure, public :: on_rejection => mh_on_rejection
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
function mh_get_chain(this, bin, err) result(rst)
    !! Gets a copy of the stored Markov chain.
    class(metropolis_hastings), intent(in) :: this
        !! The metropolis_hastings object.
    real(real64), intent(in), optional :: bin
        !! An optional input allowing for a burn-in region.  The parameter
        !! represents the amount (percentage-based) of the overall chain to 
        !! disregard as "burn-in" values.  The value shoud exist on [0, 1).
        !! The default value is 0 such that no values are disregarded.
    class(errors), intent(inout), optional, target :: err
        !! The error handling object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting chain with each parameter represented by a column.

    ! Local Variables
    integer(int32) :: npts, nvar, flag, nstart
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
    if (present(bin)) then
        nstart = floor(bin * npts)
        npts = npts - nstart
    else
        nstart = 1
    end if

    ! Process
    allocate(rst(npts, nvar), stat = flag, &
        source = this%m_buffer(nstart:,1:nvar))
    if (flag /= 0) then
        call report_memory_error(errmgr, "mh_get_chain", flag)
        return
    end if
end function

! ------------------------------------------------------------------------------
function mh_proposal(this, xc) result(rst)
    !! Proposes a new sample set of variables.  The sample is generated by
    !! sampling a multivariate normal distribution.
    class(metropolis_hastings), intent(inout) :: this
        !! The metropolis_hastings object.
    real(real64), intent(in), dimension(:) :: xc
        !! The current set of variables.
    real(real64), allocatable, dimension(:) :: rst
        !! The proposed set of variables.

    ! Sample from the distribution
    call this%m_propDist%set_means(xc) ! center the distribution at the current state
    rst = sample_normal_multivariate(this%m_propDist)
end function

! ------------------------------------------------------------------------------
function mh_hastings_ratio(this, xc, xp) result(rst)
    !! Evaluates the Hasting's ratio.  If the proposal distribution is 
    !! symmetric, this ratio is unity; however, in the case of an asymmetric
    !! distribution this ratio is not ensured to be unity.
    class(metropolis_hastings), intent(inout) :: this
        !! The metropolis_hastings object.
    real(real64), intent(in), dimension(:) :: xc
        !! The current state vector.
    real(real64), intent(in), dimension(size(xc)) :: xp
        !! The proposed state vector.
    real(real64) :: rst
        !! The ratio.

    ! Local Variables
    real(real64), allocatable, dimension(:) :: means
    real(real64) :: q1, q2

    ! Initialization
    means = this%m_propDist%get_means()

    ! Process
    call this%m_propDist%set_means(xc)
    q1 = this%m_propDist%pdf(xp)

    call this%m_propDist%set_means(xp)
    q2 = this%m_propDist%pdf(xc)

    ! Restore the state of the object
    call this%m_propDist%set_means(means)

    ! Compute the ratio
    rst = q2 / q1
end function

! ------------------------------------------------------------------------------
function mh_target(this, x) result(rst)
    !! Returns the probability value from the target distribution at the
    !! specified state.  The user is expected to overload this routine to
    !! define the desired distribution.  The default behavior of this
    !! routine is to sample a multivariate normal distribution with a mean
    !! of zero and a variance of one (identity covariance matrix).
    class(metropolis_hastings), intent(inout) :: this
        !! The metropolis_hastings object.
    real(real64), intent(in), dimension(:) :: x
        !! The state vector.
    real(real64) :: rst
        !! The value of the probability density function of the distribution
        !! being sampled.

    ! Local Variables
    integer(int32) :: i, n
    real(real64), allocatable, dimension(:) :: mu
    real(real64), allocatable, dimension(:,:) :: sigma
    type(multivariate_normal_distribution) :: dist

    ! The default will be a multivariate normal distribution with an identity
    ! matrix for the covariance matrix
    n = size(x)
    allocate(mu(n), sigma(n, n), source = 0.0d0)
    do i = 1, n
        sigma(i,i) = 1.0d0
    end do
    call dist%initialize(mu, sigma)
    rst = dist%pdf(x)
end function

! ------------------------------------------------------------------------------
subroutine mh_sample(this, xi, niter, err)
    !! Samples the distribution using the Metropolis-Hastings algorithm.
    class(metropolis_hastings), intent(inout) :: this
        !! The metropolis_hastings object.
    real(real64), intent(in), dimension(:) :: xi
        !! An N-element array containing initial starting values of the state 
        !! variables.
    integer(int32), intent(in), optional :: niter
        !! An optional input defining the number of iterations to take.  The
        !! default is 10,000.
    class(errors), intent(inout), optional, target :: err
        !! The error handling object.

    ! Local Variables
    integer(int32) :: i, n, npts, flag
    real(real64) :: r, pp, pc, a, a1, a2, alpha
    real(real64), allocatable, dimension(:) :: xc, xp, means
    real(real64), allocatable, dimension(:,:) :: sigma
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (present(niter)) then
        npts = niter
    else
        npts = this%initial_iteration_estimate
    end if
    n = size(xi)

    ! Initialize the proposal distribution.  Use an identity matrix for the
    ! covariance matrix and assume a zero mean.
    allocate(sigma(n, n), means(n), source = 0.0d0, stat = flag)
    if (flag /= 0) then
        call report_memory_error(errmgr, "mh_sample", flag)
        return
    end if
    do i = 1, n
        sigma(i,i) = 1.0d0
    end do
    call this%m_propDist%initialize(means, sigma, err = errmgr)
    if (errmgr%has_error_occurred()) return

    ! Store the initial value
    call this%push_new_state(xi, err = errmgr)
    if (errmgr%has_error_occurred()) return

    ! Iteration Process
    xc = xi
    pc = this%target_distribution(xc)
    do i = 2, npts
        ! Create a proposal & evaluate it's PDF
        xp = this%generate_proposal(xc)
        pp = this%target_distribution(xp)

        ! Evaluate the probabilities
        a1 = this%compute_hastings_ratio(xc, xp)
        a2 = pp / pc
        alpha = min(1.0d0, a1 * a2)

        ! Do we keep our current state or move to the new state?
        call random_number(r)
        if (r <= alpha) then
            ! Take the new value
            call this%push_new_state(xp, err = errmgr)
            if (errmgr%has_error_occurred()) return

            ! Update the values
            xc = xp
            pc = pp

            ! Take additional actions on success???
            call this%on_acceptance(i, alpha, xc, xp, err = errmgr)
            if (errmgr%has_error_occurred()) return
        else
            ! Keep our current estimate
            call this%push_new_state(xc, err = errmgr)
            if (errmgr%has_error_occurred()) return

            ! Take additional actions on failure???
            call this%on_rejection(i, alpha, xc, xp, err = errmgr)
            if (errmgr%has_error_occurred()) return
        end if
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine mh_clear_chain(this)
    !! Resets the object and clears out the buffer storing the chain values.
    class(metropolis_hastings), intent(inout) :: this
        !! The metropolis_hastings object.

    ! Clear the buffer
    this%m_bufferSize = 0
    this%m_numVars = 0
end subroutine

! ------------------------------------------------------------------------------
subroutine mh_on_success(this, iter, alpha, xc, xp, err)
    !! Currently, this routine does nothing and is a placeholder for the user
    !! that inherits this class to provide functionallity upon acceptance of
    !! a proposed value.
    class(metropolis_hastings), intent(inout) :: this
        !! The metropolis_hastings object.
    integer(int32), intent(in) :: iter
        !! The current iteration number.
    real(real64), intent(in) :: alpha
        !! The proposal probabilty term used for acceptance criteria.
    real(real64), intent(in), dimension(:) :: xc
        !! An N-element array containing the current state variables.
    real(real64), intent(in), dimension(size(xc)) :: xp
        !! An N-element array containing the proposed state variables that
        !! were just accepted.
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.
end subroutine

! ------------------------------------------------------------------------------
subroutine mh_on_rejection(this, iter, alpha, xc, xp, err)
    !! Currently, this routine does nothing and is a placeholder for the user
    !! that inherits this class to provide functionallity upon rejection of
    !! a proposed value.
    class(metropolis_hastings), intent(inout) :: this
        !! The metropolis_hastings object.
    integer(int32), intent(in) :: iter
        !! The current iteration number.
    real(real64), intent(in) :: alpha
        !! The proposal probabilty term used for acceptance criteria.
    real(real64), intent(in), dimension(:) :: xc
        !! An N-element array containing the current state variables.
    real(real64), intent(in), dimension(size(xc)) :: xp
        !! An N-element array containing the proposed state variables that
        !! were just rejected.
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.
end subroutine

! ------------------------------------------------------------------------------
end module
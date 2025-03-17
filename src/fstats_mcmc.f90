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
    use collections
    implicit none
    private
    public :: metropolis_hastings
    public :: mcmc_target
    public :: evaluate_model
    
    type, abstract :: mcmc_target
        !! Defines a model of the target distribution(s).
        type(list), private :: m_items
            !! The list of parameters.
        real(real64), private, allocatable, dimension(:) :: m_y
            !! A workspace array for containing the model values.
        real(real64), public :: variance_prior = 1.0d-3
            !! The target variance prior estimate.
    contains
        procedure(evaluate_model), deferred, public :: model
        procedure, public :: get_parameter_count => mt_get_param_count
        procedure, public :: add_parameter => mt_add_param
        procedure, public :: get_parameter => mt_get_param
        procedure, public :: likelihood => mt_likelihood
        procedure, public :: evaluate_variance_prior => mt_eval_var_prior
    end type

    interface
        subroutine evaluate_model(this, xdata, xc, y)
            !! Evaluates the model at the supplied values.
            use iso_fortran_env, only : real64
            import mcmc_target
            class(mcmc_target), intent(in) :: this
                !! The mcmc_target object.
            real(real64), intent(in), dimension(:) :: xdata
                !! An M-element array containing the values at which to evaluate
                !! the model.
            real(real64), intent(in), dimension(:) :: xc
                !! An N-element array containing the model parameters.
            real(real64), intent(out), dimension(:) :: y
                !! An M-element array where the resulting model values wil
                !! be written.
        end subroutine
    end interface

! ------------------------------------------------------------------------------
    type, abstract :: mcmc_proposal
        !! Defines a type responsible for generating a proposal state for a
        !! Monte-Carlo, Markov-Chain sampler.
    contains
        procedure(sample_generator), deferred, public :: generate_sample
    end type

    interface
        subroutine sample_generator(this, tgt, xc, xp, vc, vp)
            !! The signature of a subroutine meant to generate a new sample.
            use iso_fortran_env, only : real64
            import mcmc_proposal
            import mcmc_target
            class(mcmc_proposal), intent(inout) :: this
                !! The mcmc_proposal object.
            class(mcmc_target), intent(in) :: tgt
                !! The mcmc_target object.
            real(real64), intent(in), dimension(:) :: xc
                !! The current state vector.
            real(real64), intent(out), dimension(:) :: xp
                !! The proposed state vector.
            real(real64), intent(in) :: vc
                !! The current value of the variance prior.
            real(real64), intent(out) :: vp
                !! The proposed value of the variance prior.
        end subroutine
    end interface



!     type, abstract :: mcmc_target
!         !! Defines a model of the target density.
!         real(real64), private, allocatable, dimension(:) :: m_y
!             !! An array containing the evaluated model results.
!         real(real64), private, allocatable, dimension(:,:) :: m_identity
!             !! An identity matrix of appropriate size to evaluate the target
!             !! likelihood.
!         real(real64), private, allocatable, dimension(:) :: m_priors
!             !! A workspace array for the values of the priors.
!         real(real64), public :: prior_sigma = 1.0d-2
!             !! The standard deviation of the distribution representing the
!             !! variance prior.
!     contains
!         procedure(parameter_distributions), deferred, public :: priors
!         procedure(evaluate_model), deferred, public :: model
!         procedure, public :: likelihood => mt_likelihood
!         procedure, public :: variance_prior => mt_var_prior
!         procedure, public :: prior_transition_probability => mt_transition
!     end type

!     interface
!         function parameter_distributions(this, xc) result(rst)
!             !! Returns an array containing a proposed distribution for each
!             !! model parameter.
!             use fstats_distributions, only : distribution
!             use iso_fortran_env, only : real64
!             import mcmc_target
!             class(mcmc_target), intent(in) :: this
!                 !! The mcmc_target object.
!             real(real64), intent(in), dimension(:) :: xc
!                 !! An N-element array containing the current model parameters.
!             class(distribution), allocatable, dimension(:) :: rst
!                 !! The collection of distribution objects.  One for each
!                 !! model parameter.
!         end function

!         subroutine evaluate_model(this, xdata, xc, y)
!             !! Evaluates the model at the supplied values.
!             use iso_fortran_env, only : real64
!             import mcmc_target
!             class(mcmc_target), intent(in) :: this
!                 !! The mcmc_target object.
!             real(real64), intent(in), dimension(:) :: xdata
!                 !! An M-element array containing the values at which to evaluate
!                 !! the model.
!             real(real64), intent(in), dimension(:) :: xc
!                 !! An N-element array containing the model parameters.
!             real(real64), intent(out), dimension(:) :: y
!                 !! An M-element array where the resulting model values wil
!                 !! be written.
!         end subroutine
!     end interface



! ! ------------------------------------------------------------------------------
!     type, extends(mcmc_proposal) :: random_walk_proposal
!         !! Defines a random-walk type proposal generator.
!         real(real64), allocatable, dimension(:) :: setpoint
!             !! The initial setpoint for each of the parameters.
!     contains
!         procedure, public :: generate_sample => rwp_generate
!     end type








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
        type(multivariate_normal_distribution), private :: m_propDist
            !! The proposal distribution from which to draw samples.
        logical, private :: m_propDistInitialized = .false.
            !! Set to true if the proposal distribution object has been
            !! initialized; else, false.
        integer(int32), private :: m_accepted = 0
            !! The number of accepted steps.
    contains
        procedure, public :: get_state_variable_count => mh_get_nvars
        procedure, public :: get_chain_length => mh_get_chain_length
        procedure, public :: push_new_state => mh_push
        procedure, public :: get_chain => mh_get_chain
        ! procedure, public :: sample => mh_sample
        procedure, public :: reset => mh_clear_chain
        procedure, public :: on_acceptance => mh_on_success
        procedure, public :: on_rejection => mh_on_rejection
        procedure, public :: get_accepted_count => mh_get_num_accepted

        ! Private Routines
        procedure, private :: resize_buffer => mh_resize_buffer
        procedure, private :: get_buffer_length => mh_get_buffer_length
    end type

contains
! ******************************************************************************
! MCMC_TARGET
! ------------------------------------------------------------------------------
pure function mt_get_param_count(this) result(rst)
    !! Gets the number of model parameters.
    class(mcmc_target), intent(in) :: this
        !! The mcmc_target object.
    integer(int32) :: rst
        !! The parameter count.

    rst = this%m_items%count()
end function

! ------------------------------------------------------------------------------
subroutine mt_add_param(this, x, err)
    !! Adds a new model parameter.
    class(mcmc_target), intent(inout) :: this
        !! The mcmc_target object.
    class(distribution), intent(in) :: x
        !! The parameter to add.
    class(errors), intent(inout), optional, target :: err
        !! The error handler object.

    ! Process
    call this%m_items%push(x, err = err)
end subroutine

! ------------------------------------------------------------------------------
function mt_get_param(this, i) result(rst)
    !! Gets a pointer to the stored parameter.
    class(mcmc_target), intent(in) :: this
        !! The mcmc_target object.
    integer(int32), intent(in) :: i
        !! The index of the parameter to retrieve.  If outside the bounds of the
        !! collection of parameters a null pointer is returned.
    class(distribution), pointer :: rst
        !! A pointer to the requested parameter distribution.

    ! Local Variables
    class(*), pointer :: ptr

    ! Process
    ptr => this%m_items%get(i)
    rst => null()
    select type (ptr)
    class is (distribution)
        rst => ptr
    end select
end function

! ------------------------------------------------------------------------------
function mt_likelihood(this, xdata, ydata, xc, var, err) result(rst)
    !! Computes the target likelihood.
    class(mcmc_target), intent(inout) :: this
        !! The mcmc_target object.
    real(real64), intent(in), dimension(:) :: xdata
        !! An M-element array containing the independent data points.
    real(real64), intent(in), dimension(:) :: ydata
        !! An M-element array containing the dependent data points.
    real(real64), intent(in), dimension(:) :: xc
        !! An N-element array containing the model parameters.
    real(real64), intent(in) :: var
        !! An estimate of the model variance.
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.
    real(real64) :: rst
        !! The likelihood value.

    ! Local Variables
    integer(int32) :: i, m, n, flag, ep
    real(real64) :: p, temp
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    type(normal_distribution) :: dist
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    m = size(xdata)
    n = size(xc)
    dist%standard_deviation = sqrt(var)

    ! Input Checking
    if (size(ydata) /= m) then
        call report_array_size_error(errmgr, "mt_likelihood", "ydata", m, &
            size(ydata))
        return
    end if
    if (this%get_parameter_count() /= n) then
        call report_array_size_error(errmgr, "mt_likelihood", "xc", &
            this%get_parameter_count(), n)
        return
    end if

    ! Memory Allocations
    if (.not.allocated(this%m_y)) then
        allocate(this%m_y(m), stat = flag)
        if (flag /= 0) go to 10
    end if
    if (size(this%m_y) /= m) then
        deallocate(this%m_y)
        allocate(this%m_y(m), stat = flag)
        if (flag /= 0) go to 10
    end if

    ! Evaluate the model at each data point
    call this%model(xdata, xc, this%m_y)

    ! Evaluate the likelihood
    ep = 0
    temp = 1.0d0
    do i = 1, m
        ! Evaluate the distribution at y with a mean of the model result
        dist%mean_value = this%m_y(i)
        p = dist%pdf(ydata(i))
        temp = p * temp
        if (temp == 0) then
            temp = 0.0d0
            exit
        end if

        do while (abs(temp) < 1.0d0)
            temp = 1.0d1 * temp
            ep = ep - 1
        end do

        do while (abs(temp) > 1.0d1)
            temp = 1.0d-1 * temp
            ep = ep + 1
        end do
    end do
    rst = temp * (1.0d1)**temp

    ! End
    return

    ! Memory Error Handling
10  continue
    call report_memory_error(errmgr, "mt_likelihood", flag)
    return
end function

! ------------------------------------------------------------------------------
pure function mt_eval_var_prior(this, vc, x) result(rst)
    !! Evalautes the model variance prior PDF.
    class(mcmc_target), intent(in) :: this
        !! The mcmc_target object.
    real(real64), intent(in) :: vc
        !! The current value of the model variance.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the variance prior distribution PDF.
    real(real64) :: rst
        !! The value of the variance prior distribution's PDF.

    ! Local Variables
    type(log_normal_distribution) :: dist

    ! Initialization
    dist%mean_value = vc
    dist%standard_deviation = sqrt(this%variance_prior)

    ! Process
    rst = dist%pdf(x)
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------








! ! ******************************************************************************
! ! MCMC_TARGET
! ! ------------------------------------------------------------------------------
! function mt_likelihood(this, xdata, ydata, xc, var, err) result(rst)
!     !! Computes the target likelihood based upon a multivariate normal 
!     !! distribution such that \( L = N(y | f(x), \sigma^{2} I) \).
!     class(mcmc_target), intent(inout) :: this
!         !! The mcmc_target object.
!     real(real64), intent(in), dimension(:) :: xdata
!         !! An M-element array containing the independent data points.
!     real(real64), intent(in), dimension(:) :: ydata
!         !! An M-element array containing the dependent data points.
!     real(real64), intent(in), dimension(:) :: xc
!         !! An N-element array containing the model parameters.
!     real(real64), intent(in) :: var
!         !! An estimate of the model variance.
!     class(errors), intent(inout), optional, target :: err
!         !! An error handling object.
!     real(real64) :: rst
!         !! The likelihood.

!     ! Local Variables
!     integer(int32) :: i, m, flag
!     type(multivariate_normal_distribution) :: dist
!     class(errors), pointer :: errmgr
!     type(errors), target :: deferr
    
!     ! Initialization
!     if (present(err)) then
!         errmgr => err
!     else
!         errmgr => deferr
!     end if
!     m = size(xdata)

!     ! Input Checking
!     if (size(ydata) /= m) then
!         call report_array_size_error(errmgr, "mt_likelihood", "ydata", m, &
!             size(ydata))
!         return
!     end if

!     ! Local Memory Allocations - if needed
!     if (.not.allocated(this%m_y)) then
!         allocate(this%m_y(m), stat = flag)
!         if (flag /= 0) then
!             call report_memory_error(errmgr, "mt_likelihood", flag)
!             return
!         end if
!     end if

!     ! Evaluate the model
!     call this%model(xdata, xc, this%m_y)

!     ! Evaluate the normal distribution
!     if (.not.allocated(this%m_identity)) then
!         allocate(this%m_identity(m, m), stat = flag, source = 0.0d0)
!         if (flag /= 0) then
!             call report_memory_error(errmgr, "mt_likelihood", flag)
!             return
!         end if
!     end if
!     do i = 1, m
!         this%m_identity(i, i) = var
!     end do
!     call dist%initialize(this%m_y, this%m_identity, errmgr)
!     if (errmgr%has_error_occurred()) return
!     rst = dist%pdf(ydata)
! end function

! ! ------------------------------------------------------------------------------
! function mt_transition(this, xc, xp, vc, vp, err) result(rst)
!     !! Evaluates the probability of transitioning from the current set of 
!     !! parameters to the proposed set of parameters.  This is computed as
!     !! \( p(x_p | x_c) = \prod_{i=1}^{n} p(x_{p,i} | x_{c,i}) \).
!     class(mcmc_target), intent(inout) :: this
!         !! The mcmc_target object.
!     real(real64), intent(in), dimension(:) :: xc
!         !! An N-element array containing the current parameter values.
!     real(real64), intent(in), dimension(:) :: xp
!         !! An N-element array containing the proposed parameter values.
!     real(real64), intent(in) :: vc
!         !! The current value of the variance prior.
!     real(real64), intent(in) :: vp
!         !! The proposed value of the variance prior.
!     class(errors), intent(inout), optional, target :: err
!         !! An error handler object.
!     real(real64) :: rst
!         !! The resulting probability term.

!     ! Local Variables
!     real(real64) :: vpdf
!     class(distribution), allocatable :: varprior
!     class(distribution), allocatable, dimension(:) :: priors
!     class(errors), pointer :: errmgr
!     type(errors), target :: deferr
    
!     ! Initialization
!     integer(int32) :: i, n, flag, ep
!     real(real64) :: p, temp
!     if (present(err)) then
!         errmgr => err
!     else
!         errmgr => deferr
!     end if
!     n = size(xc)

!     ! Input Checking
!     if (size(xp) /= n) then
!         call report_array_size_error(errmgr, "mt_transition", "xp", n, size(xp))
!         return
!     end if

!     ! Memory Allocation
!     if (.not.allocated(this%m_priors)) then
!         allocate(this%m_priors(n), stat = flag)
!         if (flag /= 0) then
!             call report_memory_error(errmgr, "mt_transition", flag)
!             return
!         end if
!     end if

!     ! Evaluate the variance prior
!     varprior = this%variance_prior(vc)
!     vpdf = varprior%pdf(vp)

!     ! Evaluate the model parameters
!     allocate(priors(n), source = this%priors(xc))

!     ! Compute the product
!     temp = 1.0d0
!     ep = 0
!     do i = 1, n
!         p = priors(i)%pdf(xp(i))
!         temp = p * temp
!         if (temp == 0) then
!             temp = 0.0d0
!             exit
!         end if

!         do while (abs(temp) < 1.0d0)
!             temp = 1.0d1 * temp
!             ep = ep - 1
!         end do

!         do while (abs(temp) > 1.0d1)
!             temp = 1.0d-1 * temp
!             ep = ep + 1
!         end do
!     end do
!     rst = vpdf * temp * (1.0d1)**ep
! end function

! ! ------------------------------------------------------------------------------
! function mt_var_prior(this, vc) result(rst)
!     !! Evaluates the model variance prior.
!     class(mcmc_target), intent(in) :: this
!         !! The mcmc_target object.
!     real(real64), intent(in) :: vc
!         !! The current value of the model variance.
!     class(distribution), allocatable :: rst
!         !! The probability distribution of the model variance prior.

!     ! Local Variables
!     type(log_normal_distribution) :: dist

!     ! Initialization
!     dist%mean_value = vc
!     dist%standard_deviation = this%prior_sigma

!     ! Process
!     rst = dist
! end function

! ! ******************************************************************************
! ! RANDOM_WALK_PROPOSAL
! ! ------------------------------------------------------------------------------
! subroutine rwp_generate(this, tgt, xc, xp, vc, vp)
!     !! Generates a new proposal.
!     class(random_walk_proposal), intent(inout) :: this
!         !! The random_walk_proposal object.
!     class(mcmc_target), intent(in) :: tgt
!         !! The mcmc_target object.
!     real(real64), intent(in), dimension(:) :: xc
!         !! The current state vector.
!     real(real64), intent(out), dimension(:) :: xp
!         !! The proposed state vector.
!     real(real64), intent(in) :: vc
!         !! The current value of the variance prior.
!     real(real64), intent(out) :: vp
!         !! The proposed value of the variance prior.
    

!     ! Parameters
!     integer(int32), parameter :: nrand = 1

!     ! Local Variables
!     integer(int32) :: i, n
!     real(real64) :: sigma, x(nrand), rng(2), avg, mx, mn
!     class(distribution), allocatable, dimension(:) :: dist
!     class(distribution), allocatable :: vdist

!     ! Initialization
!     n = size(xc)
!     allocate(dist(n), source = tgt%priors(this%setpoint))

!     ! Cycle over each parameter
!     do i = 1, n
!         ! Sample from the proposed distribution
!         rng = dist(i)%defined_range()
!         sigma = sqrt(dist(i)%variance())
!         avg = dist(i)%mean()
!         mn = max(rng(1), avg - 3.0d0 * sigma)
!         mx = min(rng(2), avg + 3.0d0 * sigma)   ! limit the search range to a practical range
!         x = rejection_sample(dist(i), nrand, mn, mx)

!         ! Store the sampled value
!         xp(i) = xc(i) + x(1)
!     end do

!     ! Deal with the variance term
!     allocate(vdist, source = tgt%variance_prior(0.0d0))
!     sigma = sqrt(vdist%variance())
!     avg = vdist%mean()
!     rng = vdist%defined_range()
!     mn = max(rng(1), avg - 3.0d0 * sigma)
!     mx = min(rng(2), avg + 3.0d0 * sigma)
!     x = rejection_sample(vdist, nrand, mn, mx)
!     vp = vc + x(1)
! end subroutine

! ! ------------------------------------------------------------------------------

! ! ------------------------------------------------------------------------------

! ! ------------------------------------------------------------------------------























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
    integer(int32) :: npts, nvar, flag, nstart, n
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    npts = this%get_chain_length()
    n = npts
    nvar = this%get_state_variable_count()
    if (present(bin)) then
        nstart = floor(bin * npts)
        npts = npts - nstart
    else
        nstart = 1
    end if

    ! Process
    allocate(rst(npts, nvar), stat = flag, &
        source = this%m_buffer(nstart:n,1:nvar))
    if (flag /= 0) then
        call report_memory_error(errmgr, "mh_get_chain", flag)
        return
    end if
end function

! ------------------------------------------------------------------------------
! subroutine mh_sample(this, xi, niter, err)
!     !! Samples the distribution using the Metropolis-Hastings algorithm.
!     class(metropolis_hastings), intent(inout) :: this
!         !! The metropolis_hastings object.
!     real(real64), intent(in), dimension(:) :: xi
!         !! An N-element array containing initial starting values of the state 
!         !! variables.
!     integer(int32), intent(in), optional :: niter
!         !! An optional input defining the number of iterations to take.  The
!         !! default is 10,000.
!     class(errors), intent(inout), optional, target :: err
!         !! The error handling object.

!     ! Local Variables
!     integer(int32) :: i, n, npts, flag
!     real(real64) :: r, pp, pc, a, a1, a2, alpha
!     real(real64), allocatable, dimension(:) :: xc, xp, means
!     real(real64), allocatable, dimension(:,:) :: sigma
!     class(errors), pointer :: errmgr
!     type(errors), target :: deferr
    
!     ! Initialization
!     if (present(err)) then
!         errmgr => err
!     else
!         errmgr => deferr
!     end if
!     if (present(niter)) then
!         npts = niter
!     else
!         npts = this%initial_iteration_estimate
!     end if
!     n = size(xi)
!     this%m_accepted = 0

!     ! Initialize the proposal distribution.  Use an identity matrix for the
!     ! covariance matrix and assume a zero mean.
!     if (.not.this%get_proposal_initialized()) then
!         call this%initialize_proposal(n, err = errmgr)
!         if (errmgr%has_error_occurred()) return
!     end if

!     ! Store the initial value
!     call this%push_new_state(xi, err = errmgr)
!     if (errmgr%has_error_occurred()) return
!     this%m_accepted = 1

!     ! Iteration Process
!     xc = xi
!     pc = this%target_distribution(xc)
!     do i = 2, npts
!         ! Create a proposal & evaluate it's PDF
!         xp = this%generate_proposal(xc)
!         pp = this%target_distribution(xp)

!         ! Evaluate the probabilities
!         a1 = this%compute_hastings_ratio(xc, xp)
!         a2 = pp / pc
!         alpha = min(1.0d0, a1 * a2)

!         ! Do we keep our current state or move to the new state?
!         call random_number(r)
!         if (r <= alpha) then
!             ! Take the new value
!             call this%push_new_state(xp, err = errmgr)
!             if (errmgr%has_error_occurred()) return

!             ! Update the values
!             xc = xp
!             pc = pp

!             ! Log the success
!             this%m_accepted = this%m_accepted + 1

!             ! Take additional actions on success???
!             call this%on_acceptance(i, alpha, xc, xp, err = errmgr)
!             if (errmgr%has_error_occurred()) return
!         else
!             ! Keep our current estimate
!             call this%push_new_state(xc, err = errmgr)
!             if (errmgr%has_error_occurred()) return

!             ! Take additional actions on failure???
!             call this%on_rejection(i, alpha, xc, xp, err = errmgr)
!             if (errmgr%has_error_occurred()) return
!         end if
!     end do
! end subroutine

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
pure function mh_get_num_accepted(this) result(rst)
    !! Gets the number of accepted steps.
    class(metropolis_hastings), intent(in) :: this
        !! The metropolis_hastings object.
    integer(int32) :: rst
        !! The number of accepted steps.
    rst = this%m_accepted
end function

! ------------------------------------------------------------------------------
! subroutine mh_sample(this, prop, likelihood, niter, err)
! end subroutine

! ------------------------------------------------------------------------------
end module
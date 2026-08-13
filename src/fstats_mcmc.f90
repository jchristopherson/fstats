module fstats_mcmc
    use iso_fortran_env
    use fstats_distributions
    use fstats_sampling
    use fstats_errors
    use fstats_descriptive_statistics
    use fstats_types
    use fstats_regression
    use fstats_sampling
    use collections
    implicit none
    private
    public :: chain_builder
    public :: mcmc_sampler
    public :: mcmc_target
    public :: evaluate_model
    public :: mcmc_proposal
    
    type, abstract :: mcmc_target
        !! Defines a model of the target distribution(s).
        type(list), private :: m_items
            !! The list of parameters.
        real(real64), private, allocatable, dimension(:) :: m_y
            !! A workspace array for containing the model values.
        real(real64), public :: data_noise = 1.0d0
            !! A parameter representing the noise in the data.
    contains
        procedure(evaluate_model), deferred, public :: model
        procedure, public :: get_parameter_count => mt_get_param_count
        procedure, public :: add_parameter => mt_add_param
        procedure, public :: get_parameter => mt_get_param
        procedure, public :: likelihood => mt_likelihood
        procedure, public :: evaluate_variance_prior => mt_eval_var_prior
        procedure, public :: sample_variance_prior => mt_sample_var_prior
        procedure, public :: evaluate_prior => mt_eval_prior
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
    type :: mcmc_proposal
        !! Defines a type responsible for generating a proposal state for a
        !! Monte-Carlo, Markov-Chain sampler.
        logical, private :: m_recenter = .true.
            !! Allow recentering?
        real(real64), private :: m_scale = 0.1d0
            !! Global proposal scale multiplier (applied to per-parameter stddev)
        real(real64), private, allocatable, dimension(:) :: m_param_scale
            !! Optional per-parameter scale multipliers. If not allocated, m_scale is used.
    contains
        procedure, public :: generate_sample => mp_gen
        procedure, public :: get_recenter => mp_get_recenter
        procedure, public :: set_recenter => mp_set_recenter
        procedure, public :: get_scale => mp_get_scale
        procedure, public :: set_scale => mp_set_scale
        procedure, public :: get_param_scale => mp_get_param_scale
        procedure, public :: set_param_scale => mp_set_param_scale
        procedure, public :: set_param_scales => mp_set_param_scales
    end type

! ------------------------------------------------------------------------------
    type chain_builder
        !! A type allowing for the construction of chain of values.
        integer(int32), private :: initial_iteration_estimate = 10000
            !! An initial estimate at the number of allowed iterations.
        integer(int32), private :: m_bufferSize = 0
            !! The actual number of states (# of used rows) in the buffer.
        real(real64), private, allocatable, dimension(:,:) :: m_buffer
            !! The buffer where each new state is stored as a row in the matrix.
        integer(int32), private :: m_numVars = 0
            !! The number of state variables.
    contains
        procedure, public :: get_state_variable_count => cb_get_nvars
        procedure, public :: get_chain_length => cb_get_chain_length
        procedure, public :: push_new_state => cb_push
        procedure, public :: get_chain => cb_get_chain
        procedure, public :: reset => cb_clear_chain

        ! Private Routines
        procedure, private :: resize_buffer => cb_resize_buffer
        procedure, private :: get_buffer_length => cb_get_buffer_length
    end type

! ------------------------------------------------------------------------------
    type, extends(chain_builder) :: mcmc_sampler
        !! An implementation of the Metropolis-Hastings algorithm for the
        !! generation of a Markov chain.
        integer(int32), private :: m_accepted = 0
            !! The number of accepted steps.
    contains
        procedure, public :: sample => ms_sample 
        procedure, public :: on_acceptance => ms_on_success
        procedure, public :: on_rejection => ms_on_rejection
        procedure, public :: get_accepted_count => ms_get_num_accepted
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
subroutine mt_add_param(this, x)
    !! Adds a new model parameter.
    class(mcmc_target), intent(inout) :: this
        !! The mcmc_target object.
    class(distribution), intent(in) :: x
        !! The parameter to add.

    ! Process
    call this%m_items%push(x)
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
function mt_likelihood(this, xdata, ydata, xc, var) result(rst)
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
    real(real64) :: rst
        !! The likelihood value.

    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Local Variables
    integer(int32) :: i, m, n
    real(real64) :: p, temp, v
    real(real64), allocatable, dimension(:) :: resid, nrm, lognrm
    
    ! Initialization
    m = size(xdata)
    n = size(xc)

    ! Input Checking
    if (size(ydata) /= m) error stop FS_ARRAY_SIZE_ERROR
    if (this%get_parameter_count() /= n) error stop FS_ARRAY_SIZE_ERROR

    ! Ensure a positive, non-zero variance
    if (var <= 0.0d0) then
        ! Replace non-positive variance with a small epsilon to avoid log/denom issues
        v = 1.0d-12
    else
        v = var
    end if

    ! Memory Allocations
    if (.not.allocated(this%m_y)) then
        allocate(this%m_y(m))
    end if
    if (size(this%m_y) /= m) then
        deallocate(this%m_y)
        allocate(this%m_y(m))
    end if

    ! Evaluate the model at each data point
    call this%model(xdata, xc, this%m_y)

    ! Compute the likelihood assuming the residual is normally distributed
    resid = ydata - this%m_y
    ! Compute log-likelihood per point and return the log-likelihood (sum)
    lognrm = (-resid**2 / (2.0d0 * v)) - log(sqrt(2.0d0 * pi * v))
    rst = sum(lognrm)
end function

! ------------------------------------------------------------------------------
pure function mt_eval_var_prior(this, x) result(rst)
    !! Evalautes the model variance prior PDF.
    class(mcmc_target), intent(in) :: this
        !! The mcmc_target object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the variance prior distribution PDF.
    real(real64) :: rst
        !! The value of the variance prior distribution's PDF.

    ! Local Variables
    type(log_normal_distribution) :: dist
    real(real64), parameter :: neg_inf = -huge(1.0d0)
    real(real64) :: pdfv

    ! Initialization
    dist%mean_value = 0.0d0
    dist%standard_deviation = this%data_noise
    
    ! Process: return log of the variance prior PDF; handle zero/negative PDF
    pdfv = dist%pdf(x)
    if (pdfv <= 0.0d0) then
        rst = neg_inf
    else
        rst = log(pdfv)
    end if
end function

! ------------------------------------------------------------------------------
function mt_sample_var_prior(this, vc, n) result(rst)
    !! Samples the variance prior distribution for the requested number of 
    !! samples.
    class(mcmc_target), intent(inout) :: this
        !! The mcmc_target object.
    real(real64), intent(in) :: vc
        !! The prior variance term.
    integer(int32), intent(in) :: n
        !! The number of samples.
    real(real64), allocatable, dimension(:) :: rst
        !! The requested samples.

    ! Local Variables
    type(log_normal_distribution) :: dist
    real(real64) :: xmin, xmax, sigma

    ! Establish an upper bounds on the sampling region
    dist%mean_value = vc
    dist%standard_deviation = this%data_noise
    sigma = dist%standard_deviation
    if (sigma == 0.0d0) then
        sigma = 1.0d0
        dist%standard_deviation = 1.0d0
    end if
    xmax = 6.0d0 * sigma
    xmin = 1.0d-6 * xmax

    ! Process
    rst = rejection_sample(dist, n, xmin, xmax)
end function

! ------------------------------------------------------------------------------
function mt_eval_prior(this, x) result(rst)
    !! Evaluates the PDF's for each parameter and computes a probability.
    class(mcmc_target), intent(in) :: this
        !! The mcmc_target object.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the values at which to evaluate each of
        !! the N parameter PDF's.
    real(real64) :: rst
        !! The resulting probability.

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: temp, p, neg_inf
    class(distribution), pointer :: dist

    neg_inf = -huge(1.0d0)
    
    ! Initialization
    n = this%get_parameter_count()

    ! Input Check
    if (size(x) /= n) error stop FS_ARRAY_SIZE_ERROR

    ! Process - use natural-log probabilities to avoid overflow/underflow issues
    temp = 0.0d0
    do i = 1, n
        ! Evaluate the distribution
        dist => this%get_parameter(i)
        if (.not.associated(dist)) error stop FS_NULL_POINTER_ERROR
        p = dist%pdf(x(i))
        if (p <= 0.0d0) then
            ! zero-probability parameter value — log prior is -inf
            rst = neg_inf
            return
        end if
        temp = temp + log(p)
    end do
    rst = temp  ! return log-prior
end function

! ******************************************************************************
! MCMC_PROPOSAL
! ------------------------------------------------------------------------------
subroutine mp_gen(this, tgt, xc, xp, vc, vp)
    !! Creates a new sample proposal.
    class(mcmc_proposal), intent(inout) :: this
        !! The mcmc_proposal object.
    class(mcmc_target), intent(inout) :: tgt
        !! The mcmc_target object.
    real(real64), intent(in), dimension(:) :: xc
        !! An N-element array containing the existing parameter estimates.
    real(real64), intent(out), dimension(:) :: xp
        !! An N-element array where the proposed parameters will be output.
    real(real64), intent(in) :: vc
        !! The current variance (noise) term value.
    real(real64), intent(out) :: vp
        !! The proposed variance (noise) value.

    ! Parameters
    integer(int32), parameter :: nsamples = 1

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: samples(nsamples), sigma, mu, xmax, xmin, mx, mn, rng(2)
    real(real64) :: u1, u2, z, u1v, u2v, zv, pscale
    class(distribution), pointer :: dist
    
    ! Initialization
    n = tgt%get_parameter_count()

    ! Input Checking
    if (size(xp) /= n) error stop FS_ARRAY_SIZE_ERROR

    ! Use a symmetric Gaussian random-walk proposal for each parameter.
    ! Proposals are: xp(i) = xc(i) + scale * sigma_i * Normal(0,1)
    ! For the variance term we propose on log-scale: vp = vc + scale * Normal(0,1)
    do i = 1, n
        ! Get the parameter distribution for scale information
        dist => tgt%get_parameter(i)
        if (.not.associated(dist)) error stop FS_NULL_POINTER_ERROR

        ! Recenter is not needed for RW proposals; ignore recenter flag
        sigma = sqrt(max(dist%variance(), 1.0d-12))

        ! Generate a normal(0,1) using Box-Muller
        call random_number(u1)
        call random_number(u2)
        if (u1 <= 0.0d0) u1 = 1.0d-12
        z = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * acos(0.0d0) * u2)

        if (allocated(this%m_param_scale)) then
            if (i <= size(this%m_param_scale)) then
                pscale = this%m_param_scale(i)
            else
                pscale = 1.0d0
            end if
        else
            pscale = 1.0d0
        end if
        xp(i) = xc(i) + this%m_scale * pscale * sigma * z
    end do

    ! Propose log-variance (vc is the current log-variance)
    ! generate another normal variate for variance step
    call random_number(u1v)
    call random_number(u2v)
    if (u1v <= 0.0d0) u1v = 1.0d-12
    zv = sqrt(-2.0d0 * log(u1v)) * cos(2.0d0 * acos(0.0d0) * u2v)
    vp = vc + this%m_scale * zv
end subroutine

! ------------------------------------------------------------------------------
pure function mp_get_recenter(this) result(rst)
    !! Gets a value determining if the parameter distributions should be
    !! recentered about the last stored position upon sampling.
    class(mcmc_proposal), intent(in) :: this
        !! The mcmc_proposal object.
    logical :: rst
        !! True if recentering is to be allowed; else, false.

    rst = this%m_recenter
end function

! ------------------------------------------------------------------------------
subroutine mp_set_recenter(this, x)
    !! Sets a value determining if the parameter distributions should be 
    !! recentered about the last stored position upon sampling.
    class(mcmc_proposal), intent(inout) :: this
        !! The mcmc_proposal object.
    logical, intent(in) :: x
        !! True if recentering is to be allowed; else, false.

    this%m_recenter = x
end subroutine

! ------------------------------------------------------------------------------
pure function mp_get_scale(this) result(rst)
    !! Gets the proposal scale multiplier.
    class(mcmc_proposal), intent(in) :: this
    real(real64) :: rst
    rst = this%m_scale
end function

! ------------------------------------------------------------------------------
subroutine mp_set_scale(this, x)
    !! Sets the proposal scale multiplier.
    class(mcmc_proposal), intent(inout) :: this
    real(real64), intent(in) :: x
    if (x > 0.0d0) then
        this%m_scale = x
    end if
end subroutine

! ------------------------------------------------------------------------------
pure function mp_get_param_scale(this, idx) result(rst)
    class(mcmc_proposal), intent(in) :: this
    integer(int32), intent(in) :: idx
    real(real64) :: rst
    if (allocated(this%m_param_scale)) then
        if (idx >= 1 .and. idx <= size(this%m_param_scale)) then
            rst = this%m_param_scale(idx)
        else
            rst = this%m_scale
        end if
    else
        rst = this%m_scale
    end if
end function

! ------------------------------------------------------------------------------
subroutine mp_set_param_scale(this, idx, val)
    class(mcmc_proposal), intent(inout) :: this
    integer(int32), intent(in) :: idx
    real(real64), intent(in) :: val
    integer(int32) :: n
    if (val <= 0.0d0) return
    if (.not.allocated(this%m_param_scale)) then
        n = max(1_int32, idx)
        allocate(this%m_param_scale(n))
        this%m_param_scale = this%m_scale
    else
        if (idx > size(this%m_param_scale)) then
            n = idx
            this%m_param_scale = [this%m_param_scale, spread(this%m_scale, 1, n - size(this%m_param_scale))]
        end if
    end if
    this%m_param_scale(idx) = val
end subroutine

! ------------------------------------------------------------------------------
subroutine mp_set_param_scales(this, arr)
    class(mcmc_proposal), intent(inout) :: this
    real(real64), intent(in), dimension(:) :: arr
    if (.not.allocated(this%m_param_scale)) then
        allocate(this%m_param_scale(size(arr)))
    else
        if (size(this%m_param_scale) /= size(arr)) then
            deallocate(this%m_param_scale)
            allocate(this%m_param_scale(size(arr)))
        end if
    end if
    this%m_param_scale = arr
end subroutine


! ******************************************************************************
! CHAIN_BUILDER
! ------------------------------------------------------------------------------
pure function cb_get_nvars(this) result(rst)
    !! Gets the number of state variables.
    class(chain_builder), intent(in) :: this
        !! The chain_builder object.
    integer(int32) :: rst
        !! The number of state variables.

    rst = this%m_numVars
end function

! ------------------------------------------------------------------------------
pure function cb_get_chain_length(this) result(rst)
    !! Gets the length of the chain (number of stored state variables).
    class(chain_builder), intent(in) :: this
        !! The chain_builder object.
    integer(int32) :: rst
        !! The chain length.

    rst = this%m_bufferSize
end function

! ------------------------------------------------------------------------------
subroutine cb_resize_buffer(this)
    !! Resizes the buffer to accept more states.
    class(chain_builder), intent(inout) :: this
        !! The chain_builder object.

    ! Local Variables
    integer(int32) :: m, n, mOld
    real(real64), allocatable, dimension(:,:) :: copy
    
    ! Initialization
    m = this%initial_iteration_estimate
    n = this%get_state_variable_count()

    ! Is this the first time?
    if (.not.allocated(this%m_buffer)) then
        allocate(this%m_buffer(m, n))
    end if

    ! If we're here, then we need to create a copy and go from there
    m = size(this%m_buffer, 1)
    mOld = m
    allocate(copy(m, n), source = this%m_buffer)
    deallocate(this%m_buffer)
    m = m + this%initial_iteration_estimate
    allocate(this%m_buffer(m, n))
    this%m_buffer(1:mOld,:) = copy
    deallocate(copy)
end subroutine

! ------------------------------------------------------------------------------
pure function cb_get_buffer_length(this) result(rst)
    !! Gets the actual length of the buffer.  This value will likely exceed the
    !! actual number of items in the chain.
    class(chain_builder), intent(in) :: this
        !! The chain_builder object.
    integer(int32) :: rst
        !! The actual buffer length.

    if (allocated(this%m_buffer)) then
        rst = size(this%m_buffer, 1)
    else
        rst = 0
    end if
end function

! ------------------------------------------------------------------------------
subroutine cb_push(this, x)
    !! Pushes a new set of state variables onto the buffer.
    class(chain_builder), intent(inout) :: this
        !! The chain_builder object.
    real(real64), intent(in), dimension(:) :: x
        !! The new N-element state array.

    ! Local Variables
    integer(int32) :: n, n1, nbuffer, nvars
    
    ! Initialization
    n = this%get_chain_length()
    n1 = n + 1
    nbuffer = this%get_buffer_length()
    nvars = size(x)

    ! If this is the first time, ensure the collection is initialized
    if (this%get_state_variable_count() == 0) then
        this%m_numVars = nvars
    end if

    ! Input Checking
    if (nvars /= this%get_state_variable_count()) error stop FS_ARRAY_SIZE_ERROR

    ! Resize the buffer, if necessary
    if (n == 0 .or. n == nbuffer) then
        call this%resize_buffer()
    end if

    ! Store the new state
    this%m_buffer(n1,:) = x
    this%m_bufferSize = this%m_bufferSize + 1
end subroutine

! ------------------------------------------------------------------------------
pure function cb_get_chain(this, bin) result(rst)
    !! Gets a copy of the stored Markov chain.
    class(chain_builder), intent(in) :: this
        !! The chain_builder object.
    real(real64), intent(in), optional :: bin
        !! An optional input allowing for a burn-in region.  The parameter
        !! represents the amount (percentage-based) of the overall chain to 
        !! disregard as "burn-in" values.  The value shoud exist on [0, 1).
        !! The default value is 0 such that no values are disregarded.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting chain with each parameter represented by a column.

    ! Local Variables
    integer(int32) :: npts, nvar, nstart, n
    
    ! Initialization
    npts = this%get_chain_length()
    n = npts
    nvar = this%get_state_variable_count()
    if (present(bin)) then
        nstart = int(floor(bin * real(n))) + 1
        if (nstart < 1) nstart = 1
        if (nstart > n) then
            ! No samples after burn-in
            allocate(rst(0, nvar))
        end if
        npts = n - nstart + 1
    else
        nstart = 1
    end if

    ! Process
    allocate(rst(npts, nvar), source = this%m_buffer(nstart:n,1:nvar))
end function

! ------------------------------------------------------------------------------
subroutine cb_clear_chain(this)
    !! Resets the object and clears out the buffer storing the chain values.
    class(chain_builder), intent(inout) :: this
        !! The chain_builder object.

    ! Clear the buffer
    this%m_bufferSize = 0
    this%m_numVars = 0
end subroutine

! ******************************************************************************
! MCMC_SAMPLER
! ------------------------------------------------------------------------------
subroutine ms_on_success(this, iter, alpha, xc, xp)
    !! Currently, this routine does nothing and is a placeholder for the user
    !! that inherits this class to provide functionallity upon acceptance of
    !! a proposed value.
    class(mcmc_sampler), intent(inout) :: this
        !! The mcmc_sampler object.
    integer(int32), intent(in) :: iter
        !! The current iteration number.
    real(real64), intent(in) :: alpha
        !! The proposal probabilty term used for acceptance criteria.
    real(real64), intent(in), dimension(:) :: xc
        !! An N-element array containing the current state variables.
    real(real64), intent(in), dimension(size(xc)) :: xp
        !! An N-element array containing the proposed state variables that
        !! were just accepted.
end subroutine

! ------------------------------------------------------------------------------
subroutine ms_on_rejection(this, iter, alpha, xc, xp)
    !! Currently, this routine does nothing and is a placeholder for the user
    !! that inherits this class to provide functionallity upon rejection of
    !! a proposed value.
    class(mcmc_sampler), intent(inout) :: this
        !! The mcmc_sampler object.
    integer(int32), intent(in) :: iter
        !! The current iteration number.
    real(real64), intent(in) :: alpha
        !! The proposal probabilty term used for acceptance criteria.
    real(real64), intent(in), dimension(:) :: xc
        !! An N-element array containing the current state variables.
    real(real64), intent(in), dimension(size(xc)) :: xp
        !! An N-element array containing the proposed state variables that
        !! were just rejected.
end subroutine

! ------------------------------------------------------------------------------
pure function ms_get_num_accepted(this) result(rst)
    !! Gets the number of accepted steps.
    class(mcmc_sampler), intent(in) :: this
        !! The mcmc_sampler object.
    integer(int32) :: rst
        !! The number of accepted steps.
    rst = this%m_accepted
end function

! ------------------------------------------------------------------------------
subroutine ms_sample(this, xdata, ydata, prop, tgt, niter)
    !! Samples the distribution using the Metropolis-Hastings approach.
    class(mcmc_sampler), intent(inout) :: this
        !! The mcmc_sampler object.
    real(real64), intent(in), dimension(:) :: xdata
        !! An M-element array containing the independent coordinate data points.
    real(real64), intent(in), dimension(:) :: ydata
        !! An M-element array containing the dependent coordinate data points.
    class(mcmc_proposal), intent(inout) :: prop
        !! A proposal generation object.
    class(mcmc_target), intent(inout) :: tgt
        !! An mcmc_target-based object containing the model and allowing for
        !! evaluation of likelihoods.
    integer(int32), intent(in), optional :: niter
        !! An optional input defining the number of iterations to take.  The
        !! default is 10,000.

    ! Local Variables
    integer(int32) :: i, n, n1, npts, m
    integer(int32) :: last_accept_count, adapt_interval
    real(real64) :: pp, pc, alpha, r, qprior, qvar
    real(real64) :: pc_log, pp_log, delta, u, current_variance
    real(real64) :: adapt_gain, target_accept, accept_rate, cur_scale, new_scale
    real(real64), allocatable, dimension(:) :: buffer, xc
    class(distribution), pointer :: dist
    
    ! Initialization
    if (present(niter)) then
        npts = niter
    else
        npts = this%initial_iteration_estimate
    end if
    m = size(xdata)
    n = tgt%get_parameter_count()
    n1 = n + 1  ! include the variance term (we store log-variance in the last element)
    this%m_accepted = 0

    ! Input Checking
    if (size(ydata) /= m) error stop FS_ARRAY_SIZE_ERROR

    ! Memory Allocations
    allocate(buffer(n1), xc(n1), source = 0.0d0)

    ! Get an initial starting point based upon the prior means
    do i = 1, n
        dist => tgt%get_parameter(i)
        xc(i) = dist%mean()
    end do
    ! Initialize log-variance to 0 (i.e., variance = 1.0)
    xc(n1) = 0.0d0

    ! Compute initial log-posterior for the starting state
    current_variance = exp(xc(n1))
    pc_log = tgt%likelihood(xdata, ydata, xc(1:n), current_variance)
    pc_log = pc_log + tgt%evaluate_prior(xc(1:n))
    pc_log = pc_log + tgt%evaluate_variance_prior(current_variance)

    ! Store the starting location
    call this%push_new_state(xc)

    ! Setup adaptation parameters
    adapt_interval = 100
    adapt_gain = 0.05d0
    target_accept = 0.234d0
    last_accept_count = 0

    ! Process
    do i = 2, npts
        ! Create a proposal (note: xc(1:n) are parameters, xc(n1) is log-variance)
        call prop%generate_sample(tgt, xc(1:n), buffer(1:n), xc(n1), buffer(n1))

        ! Evaluate the log-likelihood for proposed state
        pp_log = tgt%likelihood(xdata, ydata, buffer(1:n), exp(buffer(n1)))

        ! Evaluate log-prior for parameters
        qprior = tgt%evaluate_prior(buffer(1:n))

        ! Evaluate log-prior for variance (input in variance-space)
        qvar = tgt%evaluate_variance_prior(exp(buffer(n1)))

        ! Sum to get proposed log-posterior
        pp_log = pp_log + qprior + qvar

        ! MH acceptance in log-domain
        delta = pp_log - pc_log
        call random_number(u)
        if (log(u) <= delta) then
            ! Accept
            call this%push_new_state(buffer)

            xc = buffer
            pc_log = pp_log

            this%m_accepted = this%m_accepted + 1

            ! Call user hook
            call this%on_acceptance(i, min(1.0d0, exp(delta)), xc, buffer)
        else
            ! Reject: keep current state
            call this%push_new_state(xc)

            call this%on_rejection(i, min(1.0d0, exp(delta)), xc, buffer)
        end if

        ! Adapt proposal scale every adapt_interval iterations (simple global adaptation)
        if (mod(i, adapt_interval) == 0) then
            accept_rate = real(this%m_accepted - last_accept_count) / real(adapt_interval)
            last_accept_count = this%m_accepted
            ! Update global scale multiplicatively toward target acceptance
            cur_scale = prop%get_scale()
            new_scale = cur_scale * exp(adapt_gain * (accept_rate - target_accept))
            call prop%set_scale(new_scale)
        end if
    end do
end subroutine

! ------------------------------------------------------------------------------
end module

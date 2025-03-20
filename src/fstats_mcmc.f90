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
    public :: mcmc_sampler
    public :: mcmc_target
    public :: evaluate_model
    public :: mcmc_proposal
    
    type, abstract :: mcmc_target
        !! Defines a model of the target distribution(s).  This type is key
        !! to the MCMC regression process.  The approach taken is to 
        !! evaluate the model provided here and evaluating its likelihood.  The
        !! likelihood is evaluated by computing the residual between the model
        !! and data, and making the assumption that the residual should be
        !! normally distributed.
        !!
        !! $$ L = \sum_{i=1}^{n} \ln{N(y_{i} | f(x_{i}, \theta), \sigma)} $$
        !!
        !! The logarithm of the results of the normal distribution are used as
        !! the scale of values can be quite extreme, especially if the model
        !! is far from the actual data; therefore, to avoid scaling induced
        !! overflow or underflow errors the logarithmic likelihood is utilized.
        !! The \(\sigma\) term results from the data_noise parameter and 
        !! is a representation of just that, the noise in the data.  The square 
        !! of this parameter can also be referred to as the variance prior.
        !! The default utilized here within is to assume the variance prior is 
        !! logarithmically distributed, and as such, is never negative-valued.
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
    contains
        procedure, public :: generate_sample => ms_gen
        procedure, public :: get_recenter => ms_get_recenter
        procedure, public :: set_recenter => ms_set_recenter
    end type

    type mcmc_sampler
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
        procedure, public :: sample => mh_sample
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

    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Local Variables
    integer(int32) :: i, m, n, flag
    real(real64) :: p, temp, v
    real(real64), allocatable, dimension(:) :: resid, nrm, lognrm
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    m = size(xdata)
    n = size(xc)

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

    ! Ensure a non-zero variance
    if (var == 0.0d0) then
        ! Warn the user that a zero variance term has been encountered
        ! call report_zero_variance_warning(errmgr, "mt_likelihood")
        v = 1.0d0
    else
        v = var
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

    ! Compute the likelihood assuming the residual is normally distributed
    resid = ydata - this%m_y
    lognrm = (-resid**2 / (2.0d0 * v)) - log(sqrt(2.0d0 * pi * v))
    rst = sum(lognrm)
    rst = exp(rst)

    ! End
    return

    ! Memory Error Handling
10  continue
    call report_memory_error(errmgr, "mt_likelihood", flag)
    return
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

    ! Initialization
    dist%mean_value = 0.0d0
    dist%standard_deviation = this%data_noise
    
    ! Process
    rst = dist%pdf(x)
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
function mt_eval_prior(this, x, err) result(rst)
    !! Evaluates the PDF's for each parameter and computes a probability.
    class(mcmc_target), intent(in) :: this
        !! The mcmc_target object.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the values at which to evaluate each of
        !! the N parameter PDF's.
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.
    real(real64) :: rst
        !! The resulting probability.

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: temp, p
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    class(distribution), pointer :: dist
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = this%get_parameter_count()

    ! Input Check
    if (size(x) /= n) then
        call report_array_size_error(errmgr, "mt_eval_prior", "x", n, size(x))
        return
    end if

    ! Process - use log prorabilities to avoid overflow/underflow issues
    temp = 0.0d0
    do i = 1, n
        ! Evaluate the distribution
        dist => this%get_parameter(i)
        p = log10(dist%pdf(x(i)))

        temp = temp + p
    end do
    rst = (1.0d1)**temp
end function

! ******************************************************************************
! MCMC_SAMPLER
! ------------------------------------------------------------------------------
subroutine ms_gen(this, tgt, xc, xp, vc, vp, err)
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
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.

    ! Parameters
    integer(int32), parameter :: nsamples = 1

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: samples(nsamples), sigma, mu, xmax, xmin, mx, mn, rng(2)
    class(distribution), pointer :: dist
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = tgt%get_parameter_count()

    ! Input Checking
    if (size(xp) /= n) then
        call report_array_size_error(errmgr, "ms_gen", "xp", n, size(xp))
        return
    end if

    ! Handle each parameter
    do i = 1, n
        ! Get the parameter distribution
        dist => tgt%get_parameter(i)
        if (.not.associated(dist)) then
            call report_null_pointer_error(errmgr, "ms_gen", "dist")
            return
        end if

        ! Recenter the distribution
        if (this%get_recenter()) call dist%recenter(xc(i))

        ! Get limit values for the parameter
        rng = dist%defined_range()
        mx = maxval(rng)
        mn = minval(rng)
        mu = dist%mean()
        sigma = sqrt(dist%variance())
        xmax = min(mu + 6.0d0 * sigma, mx)
        xmin = max(mu - 6.0d0 * sigma, mn)

        ! Sample the distribution
        samples = rejection_sample(dist, nsamples, xmin, xmax)
        xp(i) = samples(1)
    end do

    ! Handle the variance term
    samples = tgt%sample_variance_prior(vc, nsamples)
    vp = samples(1)
end subroutine

! ------------------------------------------------------------------------------
pure function ms_get_recenter(this) result(rst)
    !! Gets a value determining if the parameter distributions should be
    !! recentered about the last stored position upon sampling.
    class(mcmc_proposal), intent(in) :: this
        !! The mcmc_proposal object.
    logical :: rst
        !! True if recentering is to be allowed; else, false.

    rst = this%m_recenter
end function

! ------------------------------------------------------------------------------
subroutine ms_set_recenter(this, x)
    !! Sets a value determining if the parameter distributions should be 
    !! recentered about the last stored position upon sampling.
    class(mcmc_proposal), intent(inout) :: this
        !! The mcmc_proposal object.
    logical, intent(in) :: x
        !! True if recentering is to be allowed; else, false.

    this%m_recenter = x
end subroutine

! ******************************************************************************
! MCMC_SAMPLER
! ------------------------------------------------------------------------------
pure function mh_get_nvars(this) result(rst)
    !! Gets the number of state variables.
    class(mcmc_sampler), intent(in) :: this
        !! The mcmc_sampler object.
    integer(int32) :: rst
        !! The number of state variables.

    rst = this%m_numVars
end function

! ------------------------------------------------------------------------------
pure function mh_get_chain_length(this) result(rst)
    !! Gets the length of the chain (number of stored state variables).
    class(mcmc_sampler), intent(in) :: this
        !! The mcmc_sampler object.
    integer(int32) :: rst
        !! The chain length.

    rst = this%m_bufferSize
end function

! ------------------------------------------------------------------------------
subroutine mh_resize_buffer(this, err)
    !! Resizes the buffer to accept more states.
    class(mcmc_sampler), intent(inout) :: this
        !! The mcmc_sampler object.
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
    class(mcmc_sampler), intent(in) :: this
        !! The mcmc_sampler object.
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
    class(mcmc_sampler), intent(inout) :: this
        !! The mcmc_sampler object.
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
    class(mcmc_sampler), intent(in) :: this
        !! The mcmc_sampler object.
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
subroutine mh_clear_chain(this)
    !! Resets the object and clears out the buffer storing the chain values.
    class(mcmc_sampler), intent(inout) :: this
        !! The mcmc_sampler object.

    ! Clear the buffer
    this%m_bufferSize = 0
    this%m_numVars = 0
end subroutine

! ------------------------------------------------------------------------------
subroutine mh_on_success(this, iter, alpha, xc, xp, err)
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
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.
end subroutine

! ------------------------------------------------------------------------------
subroutine mh_on_rejection(this, iter, alpha, xc, xp, err)
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
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.
end subroutine

! ------------------------------------------------------------------------------
pure function mh_get_num_accepted(this) result(rst)
    !! Gets the number of accepted steps.
    class(mcmc_sampler), intent(in) :: this
        !! The mcmc_sampler object.
    integer(int32) :: rst
        !! The number of accepted steps.
    rst = this%m_accepted
end function

! ------------------------------------------------------------------------------
subroutine mh_sample(this, xdata, ydata, prop, tgt, niter, err)
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
    class(errors), intent(inout), optional, target :: err
        !! The error handling object.

    ! Local Variables
    integer(int32) :: i, n, n1, npts, m, flag
    real(real64) :: pp, pc, alpha, r, qprior, qvar
    real(real64), allocatable, dimension(:) :: buffer, xc
    class(distribution), pointer :: dist
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
    m = size(xdata)
    n = tgt%get_parameter_count()
    n1 = n + 1  ! include the variance term
    this%m_accepted = 0
    pc = 1.0d0

    ! Input Checking
    if (size(ydata) /= m) then
        call report_array_size_error(errmgr, "mh_sample", "ydata", m, size(ydata))
        return
    end if

    ! Memory Allocations
    allocate(buffer(n1), xc(n1), source = 0.0d0, stat = flag)
    if (flag /= 0) go to 10

    ! Get an initial starting point based upon the prior means
    do i = 1, n
        dist => tgt%get_parameter(i)
        xc(i) = dist%mean()
    end do

    ! Store the starting location
    call this%push_new_state(xc, err = errmgr)
    if (errmgr%has_error_occurred()) return

    ! Process
    do i = 2, npts
        ! Create a proposal
        call prop%generate_sample(tgt, xc(1:n), buffer(1:n), &
            xc(n1), buffer(n1), err = errmgr)
        if (errmgr%has_error_occurred()) return

        ! Evaluate the likelihood
        pp = tgt%likelihood(xdata, ydata, buffer(1:n), buffer(n1), err = errmgr)
        if (errmgr%has_error_occurred()) return

        qprior = tgt%evaluate_prior(buffer(1:n), err = errmgr)
        if (errmgr%has_error_occurred()) return

        qvar = tgt%evaluate_variance_prior(buffer(n1))

        ! Compute the MH ratio & evaluate
        pp = pp * qprior * qvar
        alpha = pp / pc
        alpha = min(1.0d0, alpha)
        call random_number(r)
        if (r <= alpha) then
            ! Keep the proposed value
            call this%push_new_state(buffer, err = errmgr)
            if (errmgr%has_error_occurred()) return

            ! Update the values
            xc = buffer
            pc = pp

            ! Log the success
            this%m_accepted = this%m_accepted + 1

            ! Anything else?
            call this%on_acceptance(i, alpha, xc, buffer, err = errmgr)
            if (errmgr%has_error_occurred()) return
        else
            ! Keep the existing and move on
            call this%push_new_state(xc, err = errmgr)
            if (errmgr%has_error_occurred()) return

            ! Anything else?
            call this%on_rejection(i, alpha, xc, buffer, err = errmgr)
            if (errmgr%has_error_occurred()) return
        end if
    end do

    ! End
    return

    ! Memory Error Handling
10  continue
    call report_memory_error(errmgr, "mh_sample", flag)
    return
end subroutine

! ------------------------------------------------------------------------------
end module
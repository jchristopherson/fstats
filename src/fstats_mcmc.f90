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
        type(multivariate_normal_distribution), private :: m_sampleDist
            !! The sampling distribution.
        real(real64), private, allocatable, dimension(:,:) :: m_jac
            !! The current Jacobian matrix estimate.
        real(real64), private :: m_targetSD = 1.0d0
            !! The standard deviation of the target normal distribution.  This
            !! value will be set based upon variances computed as part of the
            !! initialization process, or user definition during the 
            !! initialization process.
    contains
        procedure, public :: get_state_variable_count => mh_get_nvars
        procedure, public :: get_chain_length => mh_get_chain_length
        procedure, public :: push_new_state => mh_push
        procedure, public :: proposal => mh_proposal_sampler
        procedure, public :: target_density => mh_target_density
        procedure, public :: estimate_covariance => mh_create_cov_mtx
        procedure, public :: evaluate => mh_eval
        procedure, public :: get_chain => mh_get_chain
        procedure, public :: initialize => mh_init
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
subroutine mh_proposal_sampler(this, mu, x, x1, x2)
    !! A sampler capable of generating the next proposed step in the chain.  The
    !! default sampler is a multivariate Gaussian of the form.
    !!
    !! $$ f(\vec{x}) = \frac{1}{\sqrt{\left( 2 \pi \right)^n | \Sigma |}} 
    !! \exp{ \left( -\frac{1}{2} \left( \vec{x} - \vec{\mu} \right)^T 
    !! \Sigma^{-1} \left( \vec{x} - \vec{\mu} \right) \right) } $$
    !!
    !! where, \( \vec{\mu} \) is the previously accepted state and \( \vec{x} \)
    !! is the proposed state.
    class(metropolis_hastings), intent(inout) :: this
        !! The metropolis_hastings object.
    real(real64), intent(in), dimension(:) :: mu
        !! An N-element array containing the current state.
    real(real64), intent(out), dimension(size(mu)) :: x
        !! An N-element array containing the proposed state.
    real(real64), intent(in), optional, dimension(size(mu)) :: x1
        !! An optional N-element array containing the lower limits for the state 
        !! variables.  If not supplied, no limits are imposed.
    real(real64), intent(in), optional, dimension(size(mu)) :: x2
        !! An optional N-element array containing the upper limits for the state 
        !! variables.  If not supplied, no limits are imposed.

    ! Local Variables
    integer(int32) :: i, j, jmax, n
    real(real64), allocatable, dimension(:) :: u
    real(real64), allocatable, dimension(:,:) :: L
    type(normal_distribution) :: dist

    ! Initialization
    n = this%get_state_variable_count()
    allocate(u(n))

    ! Update the mean of the distribution
    call this%m_sampleDist%set_means(mu)

    ! Generate the next state by sampling the distribution
    call random_number(u)
    L = this%m_sampleDist%get_cholesky_factored_matrix()
    x = mu + matmul(L, u)

    ! Do we need to impose limits
    if (present(x1)) then
        ! lower limits
        x = max(x, x1)
    end if

    if (present(x2)) then
        ! upper limits
        x = min(x, x2)
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine mh_target_density(this, fcn, x, y, s, f, resid, dens, stop)
    !! Computes the target density given the current state.  The target density
    !! is calculated based upon a standard normal distribution whose standard
    !! deviation is specified as part of the problem initialization centered
    !! around a mean of 0.
    class(metropolis_hastings), intent(in) :: this
        !! The metropolis_hastings object.
    procedure(regression_function), pointer, intent(in) :: fcn
        !! A pointer to the function to be fitted.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the independent data points from the 
        !! experimental data.
    real(real64), intent(in), dimension(size(x)) :: y
        !! An N-element array containing the dependent data points from the
        !! experimental data.
    real(real64), intent(in), dimension(:) :: s
        !! An array containing the current state (model parameters).
    real(real64), intent(out), dimension(size(x)) :: f
        !! An N-element array where the function values at x will be written.
    real(real64), intent(out), dimension(size(x)) :: resid
        !! An N-element array where the residual will be written.
    real(real64), intent(out) :: dens
        !! The target density.
    logical, intent(out) :: stop
        !! If true, the user requested a stop to the process.

    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Local Variables
    integer(int32) :: n
    real(real64) :: se, x2
    type(normal_distribution) :: dist

    ! Initialization
    n = this%get_state_variable_count()

    ! Evaluate the function and compute the standard deviation and Chi-squared
    ! as a measure of the error of the fit.
    call fcn(x, s, f, stop)
    if (stop) return
    resid = f - y
    se = standard_deviation(resid)
    x2 = sum(resid**2 / y) ! chi-squared

    ! Compute the value of the density function
    dist%standard_deviation = this%m_targetSD
    dist%mean_value = 0.0d0
    dens = dist%pdf(x2)
end subroutine

! ------------------------------------------------------------------------------
subroutine mh_eval(this, fcn, x, y, niter, x1, x2, err)
    !! Initializes and evaluates the Metropolis-Hastings algorithm to generate 
    !! the Markov chain.
    class(metropolis_hastings), intent(inout) :: this
        !! The metropolis_hastings object.
    procedure(regression_function), pointer, intent(in) :: fcn
        !! A pointer to the function to be fitted.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the independent data points from the 
        !! experimental data.
    real(real64), intent(in), dimension(size(x)) :: y
        !! An N-element array containing the dependent data points from the
        !! experimental data.
    integer(int32), intent(in), optional :: niter
        !! An optional input that limits the number of iterations.  The default
        !! value is 10,000.
    real(real64), intent(in), optional, dimension(:) :: x1
        !! An optional M-element array containing the lower limits for the state 
        !! variables.  If not supplied, no limits are imposed.
    real(real64), intent(in), optional, dimension(:) :: x2
        !! An optional M-element array containing the upper limits for the state 
        !! variables.  If not supplied, no limits are imposed.
    class(errors), intent(inout), optional, target :: err
        !! The error handling object.

    ! Local Variables
    logical :: check
    integer(int32) :: i, n, npts, flag, maxiter
    real(real64) :: alpha, pc, pp, r, a
    real(real64), allocatable, dimension(:) :: xp, xc, f, resid
    real(real64), allocatable, dimension(:,:) :: sigma
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = this%get_state_variable_count()
    npts = size(x)
    if (present(niter)) then
        maxiter = niter
    else
        maxiter = this%initial_iteration_estimate
    end if

    ! Check to ensure the object has been initialized
    if (this%get_state_variable_count() == 0) then
    end if

    ! Memory Allocations
    allocate(xc(n), stat = flag, source = this%m_buffer(1,:))
    if (flag /= 0) go to 10
    allocate(xp(n), f(npts), resid(npts), stat = flag)
    if (flag /= 0) go to 10

    ! Evaluate the next step
    call this%proposal(xc, xp, x1, x2)

    ! Compute the target density
    call this%target_density(fcn, x, y, xc, f, resid, pc, check)

    ! Continue the iteration process now that initial values have all been
    ! established
    do i = 2, maxiter
        ! Define a new proposal
        call this%proposal(xc, xp, x1, x2)

        ! Compute the target density
        call this%target_density(fcn, x, y, xp, f, resid, pp, check)
        if (check) return

        ! Compute the acceptance ratio & see if the step is good enough
        a = pp / pc
        alpha = min(1.0d0, a)
        print *, "ALPHA = ", alpha
        print *, "PC = ", pc
        print *, "PP = ", pp
        call random_number(r)
        if (r < alpha) then
            ! Accept the step
            call this%push_new_state(xp, errmgr)
            if (errmgr%has_error_occurred()) return
            xc = xp
            pc = pp
        else
            ! Stay where we're at
            call this%push_new_state(xc, errmgr)
            if (errmgr%has_error_occurred()) return
        end if
    end do

    ! End
    return

    ! Memory Error Handler
10  continue
    call report_memory_error(errmgr, "mh_eval", flag)
    return
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
subroutine mh_init(this, fcn, x, y, mdl, tsd, err)
    use linalg, only : trace
    !! Initializes the Markov chain and sets up the metropolis_hastings object.
    class(metropolis_hastings), intent(inout) :: this
        !! The metropolis_hastings object.
    procedure(regression_function), pointer, intent(in) :: fcn
        !! A pointer to the function to fit.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the independent data points from the 
        !! experimental data.
    real(real64), intent(in), dimension(size(x)) :: y
        !! An N-element array containing the dependent data points from the
        !! experimental data.
    real(real64), intent(in), dimension(:) :: mdl
        !! An M-element array containing an initial estimate of the model 
        !! coefficients (state variables).
    real(real64), intent(in), optional :: tsd
        !! An optional input that can be used to define the target distribution
        !! standard deviation.  If not supplied, this value is set to be 5% of
        !! the square root of the largest variance term based upon the initial
        !! solution estimate.
    class(errors), intent(inout), optional, target :: err
        !! The error handling object.

    ! Local Variables
    logical :: check
    integer(int32) :: neqn, nvar, flag
    real(real64), allocatable, dimension(:,:) :: sigma
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    neqn = size(x)
    nvar = size(mdl)

    ! Initialize the Jacobian
    if (allocated(this%m_jac)) deallocate(this%m_jac)
    allocate(this%m_jac(neqn, nvar), stat = flag)
    if (flag /= 0) then
        call report_memory_error(errmgr, "mh_init", flag)
        return
    end if

    ! Create an initial covariance matrix estimate
    sigma = this%estimate_covariance(fcn, x, mdl, check)
    call this%m_sampleDist%initialize(mdl, sigma, err = errmgr)
    if (errmgr%has_error_occurred()) return

    ! Store an estimate for the target density standard deviation based upon
    ! the variance terms.
    if (present(tsd)) then
        this%m_targetSD = tsd
    else
        this%m_targetSD = 5.0d-2 * sqrt(maxval(sigma))
    end if

    ! Store the initial estimate
    call this%push_new_state(mdl, err = errmgr)
    if (errmgr%has_error_occurred()) return
end subroutine

! ------------------------------------------------------------------------------
function mh_create_cov_mtx(this, fcn, x, mdl, stop, err) result(rst)
    !! Creates an estimate of the current estimate of the covariance matrix
    !! by computing the Jacobian matrix of the system as the Jacobian matrix is
    !! effectively a linearization of the design matrix about the current state
    !! point.  The covariance matrix is then computed as \(\Sigma = 
    !! \left(J^{T} J \right)^{-1}\).
    class(metropolis_hastings), intent(inout) :: this
        !! The metropolis_hastings object.
    procedure(regression_function), pointer, intent(in) :: fcn
        !! A pointer to the regression_function to evaluate.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing x-coordinate data.
    real(real64), intent(in), dimension(:) :: mdl
        !! An array containing the current estimate of the model parameters.
    logical, intent(out) :: stop
        !! A value that the user can set in fun forcing the
        !! evaluation process to stop prior to completion.
    class(errors), intent(inout), optional, target :: err
        !! The error handling object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting covariance matrix.

    ! Local Variables
    integer(int32) :: n, nvar, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)
    nvar = size(mdl)

    ! Ensure the storage matrix for the Jacobian is initialized
    if (.not.allocated(this%m_jac)) then
        allocate(this%m_jac(n, nvar), stat = flag)
        if (flag /= 0) then
            call report_memory_error(errmgr, "mh_create_cov_mtx", flag)
            return
        end if
    end if

    ! Compute the Jacobian matrix
    call jacobian(fcn, x, mdl, this%m_jac, stop, err = errmgr)
    if (errmgr%has_error_occurred()) return
    if (stop) return

    ! Evaluate the covariance matrix
    allocate(rst(nvar, nvar), stat = flag)
    if (flag /= 0) then
        call report_memory_error(errmgr, "mh_update_cov_mtx", flag)
        return
    end if
    call covariance_matrix(this%m_jac, rst, err = errmgr)
    if (errmgr%has_error_occurred()) return
end function

! ------------------------------------------------------------------------------
! TO DO: Define burn-in limits
! end function

! TO DO: Reset chain to zero length
! TO DO: Get covariance matrix

! ------------------------------------------------------------------------------
end module
module fstats_mcmc
    use iso_fortran_env
    use fstats_distributions
    use fstats_sampling
    use ferror
    use fstats_errors
    use fstats_descriptive_statistics
    implicit none
    private
    public :: metropolis_hastings

    interface
        pure function fit_function(x, mdl) result(rst)
            use iso_fortran_env, only : real64
            !! The signature of the function to be fit.
            real(real64), intent(in), dimension(:) :: x
                !! An N-element array containing the independent variable values
                !! at which to evaluate the function.
            real(real64), intent(in), dimension(:) :: mdl
                !! The model parameters.
            real(real64), dimension(size(x)) :: rst
                !! An N-element array containing the values of the function at
                !! each value in x.
        end function
    end interface

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
        procedure(fit_function), private, pointer, nopass :: fcn
            !! The function to fit.
    contains
        procedure, public :: get_state_variable_count => mh_get_nvars
        procedure, public :: get_chain_length => mh_get_chain_length
        procedure, public :: push_new_state => mh_push
        procedure, public :: proposal => mh_proposal_sampler
        procedure, public :: target_density => mh_target_density
        procedure, public :: estimate_covariance => mh_init_covariance
        procedure, public :: evaluate => mh_eval
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
        !! The metropolis_hasings object.
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
end subroutine

! ------------------------------------------------------------------------------
subroutine mh_proposal_sampler(this, mu, x, f, x1, x2)
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
    real(real64), intent(out) :: f
        !! The value of the probability distribution function at the proposed 
        !! state centered around the current state.
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

    ! Evaluate the pdf
    f = this%m_sampleDist%pdf(x)
end subroutine

! ------------------------------------------------------------------------------
function mh_target_density(this, x, y, s) result(rst)
    !! Computes the target density given the current state.  The target density
    !! is calculated as follows.
    !!
    !! $$ f(\vec{x}) = \frac{1}{\sqrt{2 \pi}^n \sigma} \exp{\left(-\frac{1}{2} 
    !! \left(n - 1 \right) \chi^2 \right)} $$
    class(metropolis_hastings), intent(in) :: this
        !! The metropolis_hastings object.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the independent data points from the 
        !! experimental data.
    real(real64), intent(in), dimension(size(x)) :: y
        !! An N-element array containing the dependent data points from the
        !! experimental data.
    real(real64), intent(in), dimension(:) :: s
        !! An array containing the current state (model parameters).
    real(real64) :: rst
        !! The function value.

    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Local Variables
    integer(int32) :: n
    real(real64) :: se, x2
    real(real64), allocatable, dimension(:) :: f, e

    ! Evaluate the function and compute the standard deviation and Chi-squared
    ! as a measure of the error of the fit.
    n = this%get_state_variable_count()
    f = this%fcn(x, s)
    e = f - y
    se = standard_deviation(e)
    x2 = sum(e**2 / y)

    ! Compute the value of the density function
    rst = exp(-0.5d0 * (n - 1.0d0) * x2) / (sqrt(2.0d0 * pi)**n * se)
end function

! ------------------------------------------------------------------------------
function mh_init_covariance(this, x, y, mdl) result(rst)
    !! Generates an initial estimate of the covariance matrix.  The approach
    !! used by this routine is to assume the covariance matrix takes the form
    !! of an identity matrix scaled by the variance of the error of the current
    !! model state.
    class(metropolis_hastings), intent(in) :: this
        !! The metropolis_hastings object.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the independent data points from the 
        !! experimental data.
    real(real64), intent(in), dimension(size(x)) :: y
        !! An N-element array containing the dependent data points from the
        !! experimental data.
    real(real64), intent(in), dimension(:) :: mdl
        !! An initial estimate of the model coefficients.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The initial covariance matrix estimate.

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: s2
    real(real64), allocatable, dimension(:) :: f, e

    ! Evaluate the model
    f = this%fcn(x, mdl)
    e = f - y
    s2 = variance(e)

    ! Generate the matrix
    n = size(mdl)
    allocate(rst(n, n), source = 0.0d0)
    do i = 1, n
        rst(i,i) = s2
    end do
end function

! ------------------------------------------------------------------------------
! TO DO: Update covariance matrix????

! ------------------------------------------------------------------------------
subroutine mh_eval(this, fcn, x, y, mdl, niter, x1, x2, err)
    !! Initializes and evaluates the Metropolis-Hastings algorithm to generate 
    !! the Markov chain.
    class(metropolis_hastings), intent(inout) :: this
        !! The metropolis_hastings object.
    procedure(fit_function), pointer, intent(in) :: fcn
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
    integer(int32), intent(in), optional :: niter
        !! An optional input that limits the number of iterations.  The default
        !! value is 10,000.
    real(real64), intent(in), optional, dimension(size(mdl)) :: x1
        !! An optional M-element array containing the lower limits for the state 
        !! variables.  If not supplied, no limits are imposed.
    real(real64), intent(in), optional, dimension(size(mdl)) :: x2
        !! An optional M-element array containing the upper limits for the state 
        !! variables.  If not supplied, no limits are imposed.
    class(errors), intent(inout), optional, target :: err
        !! The error handling object.

    ! Local Variables
    integer(int32) :: i, n, flag, maxiter
    real(real64) :: alpha, pc, pp, r, a, a1, a2, gc, gp
    real(real64), allocatable, dimension(:) :: xp, xc
    real(real64), allocatable, dimension(:,:) :: sigma
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(mdl)
    this%fcn => fcn
    if (present(niter)) then
        maxiter = niter
    else
        maxiter = this%initial_iteration_estimate
    end if

    ! Memory Allocations
    allocate(xc(n), stat = flag, source = mdl)
    if (flag /= 0) go to 10
    allocate(xp(n), stat = flag)
    if (flag /= 0) go to 10

    ! Create an initial covariance matrix estimate & initial the proposal 
    ! distribution
    sigma = this%estimate_covariance(x, y, xc)
    call this%m_sampleDist%initialize(xc, sigma, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Store the initial estimate of the model parameters
    call this%push_new_state(mdl, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Evaluate the next step
    call this%proposal(xc, xp, gc, x1, x2)

    ! Compute the target density
    pc = this%target_density(x, y, xc)

    ! Continue the iteration process now that initial values have all been
    ! established
    do i = 2, maxiter
        ! Define a new proposal
        call this%proposal(xc, xp, gp, x1, x2)

        ! Compute the target density
        pp = this%target_density(x, y, xp)

        ! Compute the acceptance ratio & see if the step is good enough
        a1 = pp / pc
        a2 = gp / gc
        a = a1 * a2
        alpha = min(1.0d0, a)
        call random_number(r)
        if (r < alpha) then
            ! Accept the step
            call this%push_new_state(xp, errmgr)
            if (errmgr%has_error_occurred()) return
            xc = xp
            pc = pp
            gc = gp

            ! Update the covariance matrix???
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
! function mh_get_chain(this, err) result(rst)
!     !! Gets a copy of the chain

! TO DO: Define burn-in limits
! end function

! TO DO: Reset chain to zero length
! TO DO: Get covariance matrix

! ------------------------------------------------------------------------------
end module
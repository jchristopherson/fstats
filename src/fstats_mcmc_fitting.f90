module fstats_mcmc_fitting
    use iso_fortran_env
    use ieee_arithmetic
    use fstats_mcmc
    use fstats_errors
    use ferror
    use fstats_regression
    use fstats_distributions
    use fstats_sampling
    implicit none
    private
    public :: mcmc_regression

    type, extends(metropolis_hastings) :: mcmc_regression
        real(real64), public, allocatable, dimension(:) :: x
            !! The independent-variable data to fit.
        real(real64), public, allocatable, dimension(:) :: y
            !! The dependent-variable data to fit.
        procedure(regression_function), pointer, nopass :: fcn
            !! The function to fit.
        class(distribution), public, allocatable, dimension(:) :: proposals
            !! A collection of proposal generators.  There must be one 
            !! distribution object for each model parameter.  They may each
            !! be different so long as they derive from the distribution type.
        type(log_normal_distribution), public :: variance_distribution
            !! The proposal distribution representing the variance of the fit.
            !! A log-normal distribution is chosen as the values for the
            !! variance parameter are always positive by definition.
        real(real64), public, allocatable, dimension(:) :: max_parameter_values
            !! A list of upper limits for the sampling of each model parameter, 
            !! including the model variance as the last term in the array.
        real(real64), public, allocatable, dimension(:) :: min_parameter_values
            !! A list of lower limits for the sampling of each model parameter, 
            !! including the model variance as the last term in the array.

        ! -----
        ! Private Workspace Arrays
        real(real64), private, allocatable, dimension(:) :: m_f0
            !! An N-element array used for containing the current function 
            !! estimate (N = size(x)).
    contains
        procedure, public :: generate_proposal => mr_proposal
        procedure, public :: target_distribution => mr_target
        procedure, public :: covariance_matrix => mr_covariance
    end type

contains
! ------------------------------------------------------------------------------
function mr_proposal(this, xc) result(rst)
    !! Proposes a new sample set of variables.  Be sure to have defined all
    !! distributions prior to calling this routine.  If the distributions are
    !! not defined, NaN will be returned.
    class(mcmc_regression), intent(inout) :: this
        !! The mcmc_regression object.
    real(real64), intent(in), dimension(:) :: xc
        !! The current set of model parameters.  The last entry must be the
        !! variance of the fit parameter.
    real(real64), allocatable, dimension(:) :: rst
        !! The proposed set of model parameters.

    ! Parameters
    integer(int32), parameter :: nsamples = 1

    ! Local Variables
    integer(int32) :: i, n, n1
    real(real64) :: nan, sample(nsamples)

    ! Initialization
    n1 = size(xc)
    n = n1 - 1  ! -1 accounts for the variance term
    nan = ieee_value(nan, ieee_quiet_nan)
    allocate(rst(n1), source = nan)

    ! Input Checking
    if (.not.allocated(this%proposals)) return
    if (size(this%proposals) /= n) return

    ! Process
    do i = 1, n
        sample = rejection_sample(this%proposals(i), nsamples, &
            this%min_parameter_values(i), this%max_parameter_values(i))
        rst(i) = sample(1)
    end do
    sample = rejection_sample(this%variance_distribution, nsamples, &
        this%min_parameter_values(n1), this%max_parameter_values(n1))
    rst(n1) = sample(1)
end function

! ------------------------------------------------------------------------------
function mr_target(this, x) result(rst)
    !! Returns the probability value from the target distribution given the
    !! current set of model parameters.
    !!
    !! The probability value is determined as follows, assuming \(f(x)\)
    !! is the function value.
    !! $$ \prod_{i=1}^{n} p \left( y_{i} | \theta, x_{i} \right) = 
    !! \prod_{i=1}^{n} N \left(y_{i} | f(x_{i}), \sigma^2 \right) $$.
    class(mcmc_regression), intent(inout) :: this
        !! The mcmc_regression object.
    real(real64), intent(in), dimension(:) :: x
        !! The current set of model parameters, including the model variance
        !! term as the last term in the array.
    real(real64) :: rst
        !! The value of the probability density function being sampled.

    ! Local Variables
    type(normal_distribution) :: dist
    integer(int32) :: i, npts, n, n1
    real(real64) :: p
    logical :: stop

    ! Initialization
    npts = size(this%x)
    n1 = size(x)
    n = n1 - 1

    ! Evaluate the model given the current parameters
    if (.not.allocated(this%m_f0)) allocate(this%m_f0(npts))
    call this%fcn(this%x, x(1:n), this%m_f0, stop)
    if (stop) return

    ! Evaluate the probibility distribution
    rst = 0.0d0
    dist%standard_deviation = x(n1)
    do i = 1, npts
        dist%mean_value = this%m_f0(i)
        p = log(dist%pdf(this%y(i)))
        rst = p + rst ! the product becomes a sum when working with log probabilities
    end do
end function

! ------------------------------------------------------------------------------
function mr_hastings_ratio(this, xc, xp) result(rst)
    !! Evaluates the Hasting's ratio.  If the proposal distribution is 
    !! symmetric, this ratio is unity; however, in the case of an asymmetric 
    !! distribution this ratio is not ensured to be unity.
    class(mcmc_regression), intent(inout) :: this
        !! The mcmc_regression object.
    real(real64), intent(in), dimension(:) :: xc
        !! The current set of model parameters.
    real(real64), intent(in), dimension(size(xc)) :: xp
        !! The proposed set of model parameters.
    real(real64) :: rst
        !! The ratio.

    ! TO DO, for now, just return 1
    rst = 1.0d0
end function

! ------------------------------------------------------------------------------
function mr_covariance(this, xc, err) result(rst)
    !! Computes the covariance matrix for the model given the specified model
    !! parameters.
    class(mcmc_regression), intent(inout) :: this
        !! The mcmc_regression object.
    real(real64), intent(in), dimension(:) :: xc
        !! The current set of model parameters, not including the model
        !! variance term.
    class(errors), intent(inout), optional, target :: err
        !! The error handling object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The covariance matrix.

    ! Local Variables
    logical :: stop
    integer(int32) :: m, n, flag
    real(real64) :: nan
    real(real64), allocatable, dimension(:,:) :: J
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (allocated(this%x)) then
        m = size(xc)
    else
        m = 0
    end if
    n = size(xc)
    nan = ieee_value(nan, ieee_quiet_nan)
    allocate(rst(n, n), source = nan, stat = flag)
    if (flag /= 0) go to 10

    ! Quick return
    if (m == 0) return

    ! Evaluate the Jacobian matrix
    allocate(J(m, n), stat = flag)
    if (flag /= 0) go to 10
    call jacobian(this%fcn, this%x, xc, J, stop, err = errmgr)
    if (errmgr%has_error_occurred()) return
    if (stop) return

    ! Compute the covariance matrix
    call covariance_matrix(J, rst, err = errmgr)
    if (errmgr%has_error_occurred()) return
    
    ! End
    return

    ! Memory Error Handling
10  continue
    call report_memory_error(errmgr, "mr_covariance", flag)
    return
end function

! ------------------------------------------------------------------------------
end module
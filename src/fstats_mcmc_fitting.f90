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
        !! The mcmc_regression type extends the metropolis_hastings type to
        !! specifically target regression problems.  The problem is formulated
        !! such that the target distribution takes the form \(y \sim 
        !! N \left( f(x), \sigma^{2} \right) \), where \(N\) is a normal
        !! distribution with \(f(x)\) as the mean and the model variance, 
        !! \(\sigma^2\) as an additional parameter for the algorithm to find; 
        !! therefore, the parameter state vector is of length \(n + 1\), where
        !! \(n\) is the number of model parameters and the +1 accounts for the
        !! model variance term.  The model variance proposals originate from
        !! a log-normal distribution that is seperated from the multivariate
        !! normal distribution used to generate the proposals for the model
        !! parameters.  For this reason, pay attention to the specifics of 
        !! calling each routine where model parameters are required.
        real(real64), public, allocatable, dimension(:) :: x
            !! The independent-variable data to fit.
        real(real64), public, allocatable, dimension(:) :: y
            !! The dependent-variable data to fit.
        procedure(regression_function), pointer, nopass :: fcn
            !! The function to fit.
        type(log_normal_distribution), public :: variance_distribution
            !! The proposal distribution representing the variance of the fit.
            !! A log-normal distribution is chosen as the values for the
            !! variance parameter are always positive by definition.
        real(real64), public :: max_model_variance = 1.0d3
            !! The maximum value in the search region used to sample the
            !! variance_distribution object.
        real(real64), public :: min_model_variance = 1.0d-8
            !! The minimum value in the search region used to sample the
            !! variance_distribution object.

        ! -----
        ! Private Workspace Arrays
        real(real64), private, allocatable, dimension(:) :: m_f0
            !! An N-element array used for containing the current function 
            !! estimate (N = size(x)).

        ! -----
        ! Private Member Variables
        logical, private :: m_updatePropMeans = .false.
            !! True if the the proposal means should be updated based upon the 
            !! current parameter set; else, false.
    contains
        procedure, public :: generate_proposal => mr_proposal
        procedure, public :: target_distribution => mr_target
        procedure, public :: covariance_matrix => mr_covariance
        procedure, public :: get_update_proposal_means => &
            mr_get_update_prop_means
        procedure, public :: set_update_proposal_means => &
            mr_set_update_prop_means
        procedure, public :: compute_fit_statistics => mr_calc_regression_stats
    end type

contains
! ------------------------------------------------------------------------------
function mr_proposal(this, xc) result(rst)
    !! Proposes a new sample set of variables.  Be sure to have defined all
    !! distributions prior to calling this routine.  If the distributions are
    !! not defined, default implementations with unit variance will be employed.
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
    real(real64) :: samples(nsamples)

    ! Initialization
    n1 = size(xc)
    n = n1 - 1  ! -1 accounts for the variance term
    allocate(rst(n1))

    ! Update the means
    if (this%get_update_proposal_means()) then
        call this%set_proposal_means(xc(1:n))
    end if

    ! Sample the parameters
    rst(1:n) = this%metropolis_hastings%generate_proposal(xc(1:n))

    ! Sample the model variance distribution
    samples = rejection_sample(this%variance_distribution, nsamples, &
        this%min_model_variance, this%max_model_variance)
    rst(n1) = samples(1)
end function

! ------------------------------------------------------------------------------
! https://scalismo.org/docs/Tutorials/tutorial14
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
pure function mr_get_update_prop_means(this) result(rst)
    !! Gets a value indicating if the proposal means should be updated with
    !! the current parameter set.
    class(mcmc_regression), intent(in) :: this
        !! The mcmc_regression object.
    logical :: rst
        !! True if the the proposal means should be updated based upon the 
        !! current parameter set; else, false.

    rst = this%m_updatePropMeans
end function

! ------------------------------------------------------------------------------
subroutine mr_set_update_prop_means(this, x) 
    !! Sets a value indicating if the proposal means should be updated with the
    !! current parameter set.
    class(mcmc_regression), intent(inout) :: this
        !! The mcmc_regression object.
    logical, intent(in) :: x
        !! True if the the proposal means should be updated based upon the 
        !! current parameter set; else, false.

    this%m_updatePropMeans = x
end subroutine

! ------------------------------------------------------------------------------
function mr_calc_regression_stats(this, xc, alpha, err) result(rst)
    !! Calculates statistics for the quality of fit for the regression.
    class(mcmc_regression), intent(inout) :: this
        !! The mcmc_regression object.
    real(real64), intent(in), dimension(:) :: xc
        !! The model parameters.  Be sure to only include the model parameters.
        !! Do not include the model variance that is also found as part of this
        !! analysis.
    real(real64), intent(in), optional :: alpha
        !! The significance level at which to evaluate the confidence intervals.
        !! The default value is 0.05 such that a 95% confidence interval is
        !! calculated.
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.
    type(regression_statistics), allocatable, dimension(:) :: rst
        !! The resulting statistics for each parameter.

    ! Local Variables
    logical :: stop
    integer(int32) :: n, flag
    real(real64), allocatable, dimension(:) :: resid
    real(real64), allocatable, dimension(:,:) :: c
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Compute the covariance matrix
    c = this%covariance_matrix(xc, err = errmgr)
    if (errmgr%has_error_occurred()) return

    ! Evaluate the model to compute the residuals
    n = size(this%x)
    if (.not.allocated(this%m_f0)) then
        allocate(this%m_f0(n), stat = flag)
        if (flag /= 0) then
            call report_memory_error(errmgr, "mr_calc_regression_stats", flag)
            return
        end if
    end if
    call this%fcn(this%x, xc, this%m_f0, stop)
    resid = this%m_f0 - this%y

    ! Compute the statistics
    rst = calculate_regression_statistics(resid, xc, c, alpha = alpha, &
        err = errmgr)
    if (errmgr%has_error_occurred()) return
end function

! ------------------------------------------------------------------------------
end module
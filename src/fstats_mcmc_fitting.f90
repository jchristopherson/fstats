module fstats_mcmc_fitting
    use iso_fortran_env
    use ieee_arithmetic
    use fstats_mcmc
    use fstats_errors
    use ferror
    use fstats_regression
    use fstats_distributions
    use fstats_sampling
    use fstats_descriptive_statistics
    implicit none
    private
    public :: mcmc_regression

    type, extends(metropolis_hastings) :: mcmc_regression
        !! The mcmc_regression type extends the metropolis_hastings type to
        !! specifically target regression problems.  The problem is formulated
        !! such that the target distribution takes the form \(y \sim 
        !! N \left( f(x), \sigma^{2} \right) \), where \(N\) is a normal
        !! distribution with \(f(x)\) as the mean and the model variance, 
        !! \(\sigma^2\) is determined by computing the variance for the current
        !! estimate of the model.
        real(real64), public, allocatable, dimension(:) :: x
            !! The independent-variable data to fit.
        real(real64), public, allocatable, dimension(:) :: y
            !! The dependent-variable data to fit.
        procedure(regression_function), pointer, nopass :: fcn
            !! The function to fit.
        real(real64), allocatable, dimension(:) :: upper_limits
            !! An optional array that, if used, provides an upper limit to
            !! each parameter in the model.  If used, be sure this array is the
            !! same dimension as the parameter array.  If not used, leave this
            !! alone and no upper limits will be placed on the parameters.
            !! If used and the array is not sized correctly, it will be ignored.
        real(real64), allocatable, dimension(:) :: lower_limits
            !! An optional array that, if used, provides a lower limit to
            !! each parameter in the model.  If used, be sure this array is the
            !! same dimension as the parameter array.  If not used, leave this
            !! alone and no lower limits will be placed on the parameters.
            !! If used and the array is not sized correctly, it will be ignored.

        ! -----
        ! Private Workspace Arrays
        real(real64), private, allocatable, dimension(:) :: m_f0
            !! An N-element array used for containing the current function 
            !! estimate (N = size(x)).
        real(real64), private, allocatable, dimension(:) :: m_mean
            !! A NP-element array used to contain a running mean of each
            !! parameter.

        ! -----
        ! Private Member Variables
        real(real64), private :: m_dataVariance = 1.0d0
            !! The variance within the data set itself.
    contains
        procedure, public :: generate_proposal => mr_proposal
        procedure, public :: likelihood => mr_likelihood
        procedure, public :: target_distribution => mr_target
        procedure, public :: covariance_matrix => mr_covariance
        procedure, public :: compute_fit_statistics => mr_calc_regression_stats
        procedure, public :: get_data_variance => mr_get_data_variance
        procedure, public :: set_data_variance => mr_set_data_variance
        procedure, public :: on_acceptance => mr_on_success
        procedure, public :: push_new_state => mr_push
    end type

contains
! ------------------------------------------------------------------------------
function mr_proposal(this, xc) result(rst)
    !! Arguments
    class(mcmc_regression), intent(inout) :: this
        !! The mcmc_regression object.
    real(real64), intent(in), dimension(:) :: xc
        !! The current model parameters.
    real(real64), allocatable, dimension(:) :: rst
        !! The proposed set of model parameters.

    ! Local Variables
    integer(int32) :: n

    ! Establish the parameter guess
    rst = this%metropolis_hastings%generate_proposal(xc)

    ! Apply limits?
    n = size(xc)
    if (allocated(this%upper_limits)) then
        if (size(this%upper_limits) == n) then
            rst = min(this%upper_limits, rst)
        end if
    end if
    if (allocated(this%lower_limits)) then
        if (size(this%lower_limits) == n) then
            rst = max(this%lower_limits, rst)
        end if
    end if
end function

! ------------------------------------------------------------------------------
function mr_likelihood(this, x) result(rst)
    !! Estimates the likelihood of the model.
    !!
    !! The likelihood is computed as follows assuming \(\sigma^2\) is known
    !! a priori.
    !! $$ L \left( \theta \right) = \prod_{i=1}^{n} N \left(y_{i} | f(x_{i}), 
    !! \sigma^2 \right) $$ 
    class(mcmc_regression), intent(inout) :: this
        !! The mcmc_regression object.
    real(real64), intent(in), dimension(:) :: x
        !! The current set of model parameters.
    real(real64) :: rst
        !! The likelihood value.

    ! Local Variables
    type(normal_distribution) :: dist
    integer(int32) :: i, npts, ep
    real(real64) :: p, temp
    logical :: stop

    ! Initialization
    npts = size(this%x)

    ! Evaluate the model given the current parameters
    if (.not.allocated(this%m_f0)) allocate(this%m_f0(npts))
    call this%fcn(this%x, x, this%m_f0, stop)
    if (stop) return

    ! Evaluate the probibility distribution
    temp = 1.0d0
    ep = 0
    rst = 1.0d0
    dist%standard_deviation = sqrt(this%get_data_variance())
    do i = 1, npts
        dist%mean_value = this%m_f0(i)
        p = dist%pdf(this%y(i))
        temp = p * temp
        if (temp == 0.0d0) then
            rst = 0.0d0
            return
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
    rst = temp * (1.0d1)**ep
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
    !! p \left( \theta, \sigma^2 \right) 
    !! \prod_{i=1}^{n} N \left(y_{i} | f(x_{i}), \sigma^2 \right) $$.
    class(mcmc_regression), intent(inout) :: this
        !! The mcmc_regression object.
    real(real64), intent(in), dimension(:) :: x
        !! The current set of model parameters.
    real(real64) :: rst
        !! The value of the probability density function being sampled.

    ! Process
    rst = this%likelihood(x) * this%evaluate_proposal_pdf(x)
end function

! ------------------------------------------------------------------------------
function mr_covariance(this, xc, err) result(rst)
    !! Computes the covariance matrix for the model given the specified model
    !! parameters.
    class(mcmc_regression), intent(inout) :: this
        !! The mcmc_regression object.
    real(real64), intent(in), dimension(:) :: xc
        !! The current set of model parameters.
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
        m = size(this%x)
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
function mr_calc_regression_stats(this, xc, alpha, err) result(rst)
    !! Calculates statistics for the quality of fit for the regression.
    class(mcmc_regression), intent(inout) :: this
        !! The mcmc_regression object.
    real(real64), intent(in), dimension(:) :: xc
        !! The model parameters.  Be sure to only include the model parameters.
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
pure function mr_get_data_variance(this) result(rst)
    !! Gets the variance of the observed data.
    class(mcmc_regression), intent(in) :: this
        !! The mcmc_regression object.
    real(real64) :: rst
        !! The variance.

    rst = this%m_dataVariance
end function

! ------------------------------------------------------------------------------
subroutine mr_set_data_variance(this, x)
    !! Sets the variance of the observed data.
    class(mcmc_regression), intent(inout) :: this
        !! The mcmc_regression object.
    real(real64), intent(in) :: x
        !! The variance.

    this%m_dataVariance = x
end subroutine

! ------------------------------------------------------------------------------
subroutine mr_on_success(this, iter, alpha, xc, xp, err)
    !! Updates the covariance matrix of the proposal distribution upon a 
    !! successful step.  If overloaded, be sure to call the base method to
    !! retain the functionallity required to keep the covariance matrix 
    !! up-to-date.
    class(mcmc_regression), intent(inout) :: this
        !! The mcmc_regression object.
    integer(int32), intent(in) :: iter
        !! The current iteration number.
    real(real64), intent(in) :: alpha
        !! The proposal probability term used for acceptance criteria.
    real(real64), intent(in), dimension(:) :: xc
        !! The current model parameter estimates.
    real(real64), intent(in), dimension(size(xc)) :: xp
        !! The recently accepted model parameter estimates.
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.

    ! Local Variables
    integer(int32) :: i, j, n, np
    real(real64) :: nm1, nm2, ratio
    real(real64), allocatable, dimension(:) :: delta
    real(real64), allocatable, dimension(:,:) :: sig

    ! Updates the estimate of the covariance matrix by implementing Roberts &
    ! Rosenthals adaptive approach.
    !
    ! Parameters:
    ! - xp: NP-by-1 array of the newest sampled points
    ! - xm: NP-by-1 array of the updated mean over all samples
    ! - sig: NP-by-NP old covariance matrix
    ! - n: # of samples drawn
    !
    ! C = (n - 2) / (n - 1) * sig + matmul(xp - xm, transpose(xp - xm)) / (n - 1)
    np = size(xc)
    n = this%get_chain_length()
    if (n == 1 .or. .not.allocated(this%m_mean)) then
        ! No action is necessary
        return
    end if
    nm1 = n - 1.0d0
    nm2 = n - 2.0d0
    ratio = nm2 / nm1
    delta = xp - this%m_mean
    sig = this%get_proposal_covariance()

    do j = 1, np
        do i = 1, np
            sig(i,j) = ratio * sig(i,j) + delta(i) * delta(j) / nm1
        end do
    end do

    ! Update the covariance matrix
    call this%set_proposal_covariance(sig, err = err)
end subroutine

! ------------------------------------------------------------------------------
subroutine mr_push(this, x, err)
    !! Pushes a new set of parameters onto the end of the chain buffer.
    class(mcmc_regression), intent(inout) :: this
        !! The mcmc_regression object.
    real(real64), intent(in), dimension(:) :: x
        !! The new N-element state array.
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.

    ! Local Variables
    integer(int32) :: n, npts, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Push the item onto the stack using the base method
    call this%metropolis_hastings%push_new_state(x, err = errmgr)
    if (errmgr%has_error_occurred()) return

    ! Update the running average term
    n = size(x)
    npts = this%get_chain_length()
    if (.not.allocated(this%m_mean)) then
        allocate(this%m_mean(n), stat = flag, source = x)
        if (flag /= 0) then
            call report_memory_error(errmgr, "mr_push", flag)
            return
        end if
        
        ! No more action is necessary - end here
        return
    end if

    ! Update the mean
    this%m_mean = (npts * this%m_mean + x) / (npts + 1.0d0)
end subroutine

! ------------------------------------------------------------------------------
end module
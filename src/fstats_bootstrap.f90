module fstats_bootstrap
    use iso_fortran_env
    use fstats_errors
    use omp_lib
    use fstats_distributions
    use fstats_descriptive_statistics
    use fstats_special_functions
    use fstats_regression
    implicit none
    private
    public :: bootstrap_regression_statistics
    public :: bootstrap_linear_least_squares
    public :: bootstrap_nonlinear_least_squares

! REFERENCES:
! - https://medium.com/@m21413108/bootstrapping-maximum-entropy-non-parametric-boot-python-3b1e23ea589d
! - https://cran.r-project.org/web/packages/meboot/vignettes/meboot.pdf
! - https://gist.github.com/christianjauregui/314456688a3c2fead43a48be3a47dad6

    type bootstrap_regression_statistics
        !! A container for regression-related statistical information as 
        !! computed in a bootstrap, or equivalent, calculation.
        real(real64) :: standard_error
            !! The standard error for the model coefficient.
        real(real64) :: t_statistic
            !! The T-statistic for the model coefficient.
            !!
            !! $$ t_o = \frac{ \beta_{i} }{E_{s}(\beta_{i})} $$
        real(real64) :: probability
            !! The probability that the coefficient is not statistically 
            !! important.  A statistically important coefficient will have a 
            !! low probability (p-value), typically 0.05 or lower; however, a 
            !! p-value of up to ~0.2 may be acceptable dependent upon the 
            !! problem.  Typically any p-value larger than ~0.2 indicates the 
            !! parameter is not statistically important for the model.
            !!
            !! $$ p = t_{|t_o|, df_{residual}} $$
        real(real64) :: upper_confidence_interval
            !! The upper limit of the confidence interval for the parameter.
        real(real64) :: lower_confidence_interval
            !! The lower limit of the confidence interval for the parameter.
    end type

contains
! ******************************************************************************
! LINEAR REGRESSION
! ------------------------------------------------------------------------------
subroutine bootstrap_linear_least_squares(order, intercept, x, y, &
    coeffs, ymod, resid, nsamples, stats, bias, alpha, err)
    !! Computes a linear least-squares regression to fit a set of data.
    !! Bootstrapping is utilized to gain insight into the quality of 
    !! the fit.
    integer(int32), intent(in) :: order
        !! The order of the equation to fit.  This value must be at 
        !! least one (linear equation), but can be higher as desired, 
        !! as long as there is sufficient data.
    logical, intent(in) :: intercept
        !! Set to true if the intercept is being computed as part of 
        !! the regression; else, false.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the independent variable
        !! measurement points.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the dependent variable
        !! measurement points.
    real(real64), intent(out), dimension(:) :: coeffs
        !! An ORDER+1 element array where the coefficients will
        !! be written.
    real(real64), intent(out), dimension(:) :: ymod
        !! An N-element array where the modeled data will be written.
    real(real64), intent(out), dimension(:) :: resid
        !! An N-element array where the residual error data will be 
        !! written (modeled - actual).
    integer(int32), intent(in), optional :: nsamples
        !! The number of bootstrapping samples to utilize.  
    type(bootstrap_regression_statistics), intent(out), optional, &
        dimension(:) :: stats
        !! An M-element array of bootstrap_regression_statistics items 
        !! where M = ORDER + 1 when intercept is set to true; however, 
        !! if intercept is set to false, M = ORDER.
    real(real64), intent(out), optional, dimension(:) :: bias
        !! An ORDER+1 element array where an estimate of the bias of
        !! each coefficient is returned based upon the results of the
        !! bootstrapping analysis.
    real(real64), intent(in), optional :: alpha
        !! The significance level at which to evaluate the confidence 
        !! intervals.  The default value is 0.05 such that a 95% 
        !! confidence interval is calculated.
    class(errors), intent(inout), optional, target :: err
        !! A mechanism for communicating errors and warnings to the 
        !! caller.  Possible warning and error codes are as follows.
        !! - FS_NO_ERROR: No errors encountered.
        !! - FS_ARRAY_SIZE_ERROR: Occurs if any of the arrays are not 
        !!      approriately sized.
        !! - FS_INVALID_INPUT_ERROR: Occurs if order is less than 1.
        !! - FS_MEMORY_ERROR: Occurs if there is a memory allocation 
        !!      error.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: p05 = 5.0d-2
    real(real64), parameter :: half = 5.0d-1

    ! Local Variables
    integer(int32) :: i, j, n, ns, nc, ncoeffs, flag
    real(real64) :: eps, alph
    real(real64), allocatable, dimension(:) :: fLocal, yLocal, rLocal
    real(real64), allocatable, dimension(:,:) :: allcoeffs
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (present(nsamples)) then
        ns = nsamples
    else
        ns = 1000
    end if
    if (present(alpha)) then
        alph = alpha
    else
        alph = p05
    end if
    n = size(x)
    ncoeffs = order + 1
    nc = order
    if (intercept) nc = nc + 1

    ! Compute the fit
    call linear_least_squares(order, intercept, x, y, coeffs, &
        ymod, resid, alpha = alph, err = errmgr)
    if (errmgr%has_error_occurred()) return
    
    ! Memory Allocations
    allocate(allcoeffs(ncoeffs, ns), source = zero, stat = flag)
    if (flag /= 0) then
        call report_memory_error(errmgr, "bootstrap_linear_least_squares", flag)
        return
    end if
    allcoeffs(:,1) = coeffs

    ! Establish the epsilon term
    eps = standard_deviation(y) / sqrt(real(n))
    call random_init(.false., .true.)

    ! Cycle over each data set and perform the fit
#ifdef USEOPENMP
!$OMP PARALLEL DO PRIVATE(fLocal, yLocal, rLocal) SHARED(allcoeffs)
    do i = 2, ns
        ! Allocate local arrays on a per-thread basis
        if (.not.allocated(fLocal)) allocate(fLocal(n))
        if (.not.allocated(yLocal)) allocate(yLocal(n))
        if (.not.allocated(rLocal)) allocate(rLocal(n))

        ! Compute a random data set
        call random_number(yLocal)
        yLocal = eps * (yLocal - half) + y

        ! Compute the fit of the perturbed data set
        call linear_least_squares(order, intercept, x, yLocal, &
            allcoeffs(:,i), fLocal, rLocal, alpha = alph)
    end do
!$OMP END PARALLEL DO
#else
    ! OpenMP is not available - run in a serial manner
    allocate(fLocal(n), yLocal(n), rLocal(n))
    do i = 2, ns
        ! Compute a random data set
        call random_number(yLocal)
        yLocal = eps * (yLocal - half) + y

        ! Compute the fit of the perturbed data set
        call linear_least_squares(order, intercept, x, yLocal, &
            allcoeffs(:,i), fLocal, rLocal, alpha = alph)
    end do
#endif
    
    ! Perform statistics calculations, if needed
    if (present(stats)) then
        call compute_stats(coeffs, allcoeffs, alph, intercept, stats)
    end if

    ! Compute the bias for each parameter, if needed
    if (present(bias)) then
        ! Verify the size of the array
        if (size(bias) /= ncoeffs) then
            call report_array_size_error(errmgr, &
                "bootstrap_linear_least_squares", "bias", ncoeffs, size(bias))
            return
        end if

        ! Perform the calculations
        do i = 1, ncoeffs
            bias(i) = coeffs(i) - mean(allcoeffs(i,:))
        end do
    end if
end subroutine

! ******************************************************************************
! NONLINEAR REGRESSION
! ------------------------------------------------------------------------------
subroutine bootstrap_nonlinear_least_squares(fun, x, y, params, ymod, resid, &
    nsamples, weights, maxp, minp, stats, alpha, controls, settings, info, &
    bias, err)
    !! Performs a nonlinear regression to fit a model using a version
    !! of the Levenberg-Marquardt algorithm.  Bootstrapping is utilized to gain 
    !! insight into the quality of the fit.
    procedure(regression_function), intent(in), pointer :: fun
        !! A pointer to the regression_function to evaluate.
    real(real64), intent(in) :: x(:)
        !! The M-element array containing independent data.
    real(real64), intent(in) :: y(:)
        !! The M-element array containing dependent data.
    real(real64), intent(inout) :: params(:)
        !! On input, the N-element array containing the initial estimate
        !! of the model parameters.  On output, the computed model 
        !! parameters.
    real(real64), intent(out) :: ymod(:)
        !! An M-element array where the modeled dependent data will
        !! be written.
    real(real64), intent(out) :: resid(:)
        !! An M-element array where the model residuals will be
        !! written.
    integer(int32), intent(in), optional :: nsamples
        !! The number of bootstrapping samples to utilize.  
    real(real64), intent(in), optional, target :: weights(:)
        !! An optional M-element array allowing the weighting of
        !! individual points.
    real(real64), intent(in), optional, target :: maxp(:)
        !! An optional N-element array that can be used as upper limits 
        !! on the parameter values.  If no upper limit is requested for
        !! a particular parameter, utilize a very large value.  The 
        !! internal default is to utilize huge() as a value.
    real(real64), intent(in), optional, target :: minp(:)
        !! An optional N-element array that can be used as lower limits 
        !! on the parameter values.  If no lower limit is requested for
        !! a particalar parameter, utilize a very large magnitude, but 
        !! negative, value.  The internal default is to utilize -huge() 
        !! as a value.
    type(bootstrap_regression_statistics), intent(out), optional :: stats(:)
        !! An optional N-element array that, if supplied, will be used 
        !! to return statistics about the fit for each parameter.
    real(real64), intent(in), optional :: alpha
        !! The significance level at which to evaluate the confidence 
        !! intervals.  The default value is 0.05 such that a 95% 
        !! confidence interval is calculated.
    type(iteration_controls), intent(in), optional :: controls
        !! An optional input providing custom iteration controls.
    type(lm_solver_options), intent(in), optional :: settings
        !! An optional input providing custom settings for the solver.
    type(convergence_info), intent(out), optional, target :: info
        !! An optional output that can be used to gain information about
        !! the iterative solution and the nature of the convergence.
    real(real64), intent(out), optional, dimension(:) :: bias
        !! An optional N-element array that, if supplied, will be used to 
        !! provide an estimate of the bias of each model parameter based upon
        !! the results of the bootstrapping analysis.
    class(errors), intent(inout), optional, target :: err
        !! A mechanism for communicating errors and warnings to the 
        !! caller.  Possible warning and error codes are as follows.
        !! - FS_NO_ERROR: No errors encountered.
        !! - FS_ARRAY_SIZE_ERROR: Occurs if any of the arrays are not 
        !!      properly sized.
        !! - FS_MEMORY_ERROR: Occurs if there is a memory allocation 
        !!      error.
        !! - FS_UNDERDEFINED_PROBLEM_ERROR: Occurs if the problem posed 
        !!      is underdetetermined (M < N).
        !! - FS_TOLERANCE_TOO_SMALL_ERROR: Occurs if any supplied 
        !!      tolerances are too small to be practical.
        !! - FS_TOO_FEW_ITERATION_ERROR: Occurs if too few iterations 
        !!      are allowed.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: p05 = 5.0d-2
    real(real64), parameter :: half = 5.0d-1
    
    ! Local Variables
    integer(int32) :: i, n, ns, nparams, flag
    real(real64) :: eps, alph
    real(real64), allocatable, dimension(:) :: fLocal, yLocal, rLocal
    real(real64), allocatable, dimension(:,:) :: allcoeffs
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (present(nsamples)) then
        ns = nsamples
    else
        ns = 1000
    end if
    if (present(alpha)) then
        alph = alpha
    else
        alph = p05
    end if
    n = size(x)
    nparams = size(params)

    ! Compute the fit
    call nonlinear_least_squares(fun, x, y, params, ymod, resid, &
        weights = weights, maxp = maxp, minp = minp, alpha = alph, &
        controls = controls, settings = settings, info = info, err = err)

    ! Memory Allocations
    allocate(allcoeffs(nparams, ns), source = zero, stat = flag)
    if (flag /= 0) then
        call report_memory_error(errmgr, "bootstrap_nonlinear_least_squares", &
            flag)
        return
    end if
    allcoeffs(:,1) = params
    
    ! Define initial guesses for each step.  Base upon the results of the 
    ! initial analysis as this should provide a strong starting point for
    ! subsequent analysis
    do i = 1, nparams
        allcoeffs(i,:) = params(i)
    end do

    ! Establish the epsilon term
    eps = standard_deviation(y) / sqrt(real(n))
    call random_init(.false., .true.)

    ! Cycle over each data set and perform the fit
#ifdef USEOPENMP
!$OMP PARALLEL DO
    do i = 2, ns
        ! Allocate local arrays on a per-thread basis
        if (.not.allocated(fLocal)) allocate(fLocal(n))
        if (.not.allocated(yLocal)) allocate(yLocal(n))
        if (.not.allocated(rLocal)) allocate(rLocal(n))

        ! Compute a random data set
        call random_number(yLocal)
        yLocal = eps * (yLocal - half) + y

        ! Compute the fit of the perturbed data set
        call nonlinear_least_squares(fun, x, yLocal, allcoeffs(:,i), fLocal, &
            rLocal, weights = weights, maxp = maxp, minp = minp, alpha = alph, &
            controls = controls, settings = settings, info = info)
    end do
!$OMP END PARALLEL DO
#else
    ! OpenMP is not available - run in a serial manner
    allocate(fLocal(n), yLocal(n), rLocal(n))
    do i = 2, ns
        ! Compute a random data set
        call random_number(yLocal)
        yLocal = eps * (yLocal - half) + y

        ! Compute the fit of the perturbed data set
        call nonlinear_least_squares(fun, x, yLocal, allcoeffs(:,i), fLocal, &
            rLocal, weights = weights, maxp = maxp, minp = minp, alpha = alph, &
            controls = controls, settings = settings, info = info)
    end do
#endif

    ! Perform the statistics calculations, if needed
    if (present(stats)) then
        call compute_stats(params, allcoeffs, alph, .true., stats)
    end if

    ! Compute the bias for each parameter, if needed
    if (present(bias)) then
        ! Verify the size of the array
        if (size(bias) /= nparams) then
            call report_array_size_error(errmgr, &
                "bootstrap_nonlinear_least_squares", "bias", &
                nparams, size(bias))
            return
        end if

        ! Perform the calculations
        do i = 1, nparams
            bias(i) = params(i) - mean(allcoeffs(i,:))
        end do
    end if
end subroutine

! ******************************************************************************
! PRIVATE HELPER ROUTINES
! ------------------------------------------------------------------------------
subroutine compute_stats(mdl, coeffs, alpha, intercept, stats)
    ! Arguments
    real(real64), intent(in), dimension(:) :: mdl
    real(real64), intent(inout), dimension(:,:) :: coeffs
    real(real64), intent(in) :: alpha
    logical, intent(in) :: intercept
    type(bootstrap_regression_statistics), intent(out), dimension(:) :: stats

    ! Parameters
    real(real64), parameter :: half = 0.5d0

    ! Local Variables
    integer(int32) :: i, j, i1, i2, ncoeffs, nc, nsamples
    real(real64) :: ms
    type(t_distribution) :: dist

    ! Initialization
    ncoeffs = size(coeffs, 1)
    nsamples = size(coeffs, 2)
    nc = ncoeffs
    if (.not.intercept) nc = ncoeffs - 1
    i1 = floor(half * alpha * nsamples, int32)
    i2 = nsamples - i1 + 1
    dist%dof = real(nsamples - nc)

    ! Process
    j = 1
    if (intercept) j = 0
    do i = 1, nc
        j = j + 1
        ms = trimmed_mean(coeffs(j,:), p = half * alpha)

        ! As we have a distribution of mean values, the standard deviation
        ! of this population yields the standard error estimate for the
        ! overall problem
        stats(i)%standard_error = standard_deviation(coeffs(j,:))

        ! As before, this is a distribution of mean values.  The CI can
        ! be directly estimated by considering the values of the bottom
        ! alpha/2 and top alpha/2 terms.
        stats(i)%upper_confidence_interval = coeffs(j,i2)
        stats(i)%lower_confidence_interval = coeffs(j,i1)

        ! Compute the remaining parameters
        stats(i)%t_statistic = mdl(j) / stats(i)%standard_error
        stats(i)%probability = regularized_beta(half * dist%dof, half, &
            dist%dof / (dist%dof + (stats(i)%t_statistic)**2))
    end do
end subroutine

! ------------------------------------------------------------------------------
end module
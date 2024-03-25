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
    public :: bootstrap_resampling_routine
    public :: bootstrap_statistic_routine
    public :: random_resample
    public :: bootstrap_statistics
    public :: bootstrap
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

    type bootstrap_statistics
        !! A collection of statistics resulting from the bootstrap process.
        real(real64) :: statistic_value
            !! The value of the statistic of interest.
        real(real64) :: upper_confidence_interval
            !! The upper confidence limit on the statistic.
        real(real64) :: lower_confidence_interval
            !! The lower confidence limit on the statistic.
        real(real64) :: bias
            !! The bias in the statistic.
        real(real64) :: standard_error
            !! The standard error of the statistic.
        real(real64), allocatable, dimension(:) :: population
            !! An array of the population values generated by the bootstrap
            !! process.
    end type

    interface
        subroutine bootstrap_resampling_routine(x, xn)
            !! Defines the signature of a subroutine used to compute a 
            !! resampling of data for bootstrapping purposes.
            use iso_fortran_env, only : real64
            real(real64), intent(in), dimension(:) :: x
                !! The N-element array to resample.
            real(real64), intent(out), dimension(size(x)) :: xn
                !! An N-element array where the resampled data set will be 
                !! written.
        end subroutine

        function bootstrap_statistic_routine(x) result(rst)
            !! Defines the signature of a function for computing the desired
            !! bootstrap statistic.
            use iso_fortran_env, only : real64
            real(real64), intent(in), dimension(:) :: x
                !! The array of data to analyze.
            real(real64) :: rst
                !! The resulting statistic.
        end function
    end interface

contains
! ******************************************************************************
! RESAMPLING
! ------------------------------------------------------------------------------
subroutine random_resample(x, xn)
    !! Random resampling, with replacement, based upon a normal distribution.
    real(real64), intent(in), dimension(:) :: x
        !! The N-element array to resample.
    real(real64), intent(out), dimension(size(x)) :: xn
        !! An N-element array where the resampled data set will be written.

    ! Parameters
    real(real64), parameter :: half = 0.5d0

    ! Local Variables
    integer(int32) :: n
    real(real64) :: eps

    ! Process
    n = size(x)
    eps = standard_deviation(x) / sqrt(real(n, real64))
    call random_number(xn)
    xn = eps * (xn - half) + x
end subroutine

! ******************************************************************************
! BOOTSTRAPPING
! ------------------------------------------------------------------------------
function bootstrap(stat, x, method, nsamples, alpha) result(rst)
    !! Performs a bootstrap calculation on the supplied data set for the given
    !! statistic.  The default implementation utlizes a random resampling with 
    !! replacement.  Other resampling methods may be defined by specifying an 
    !! appropriate routine by means of the method input.
    procedure(bootstrap_statistic_routine), pointer, intent(in) :: stat
        !! The routine used to compute the desired statistic.
    real(real64), intent(in), dimension(:) :: x
        !! The N-element data set.
    procedure(bootstrap_resampling_routine), pointer, intent(in), optional :: method
        !! An optional pointer to the method to use for resampling of the data.
        !! If no method is supplied, a random resampling is utilized.
    integer(int32), intent(in), optional :: nsamples
        !! An optional input, that if supplied, specifies the number of 
        !! resampling runs to perform.  The default is 10 000.
    real(real64), intent(in), optional :: alpha
        !! An optional input, that if supplied, defines the significance level
        !! to use for the analysis.  The default is 0.05.
    type(bootstrap_statistics) :: rst
        !! The resulting bootstrap_statistics type containing the confidence
        !! intervals, bias, standard error, etc. for the analyzed statistic.

    ! Parameters
    real(real64), parameter :: half = 0.5d0
    real(real64), parameter :: p05 = 5.0d-2

    ! Local Variables
    integer(int32) :: i, i1, i2, n, ns
    real(real64) :: a
    real(real64), allocatable, dimension(:) :: xn
    procedure(bootstrap_resampling_routine), pointer :: resample

    ! Initialization
    n = size(x)
    if (present(method)) then
        resample => method
    else
        resample => random_resample
    end if
    if (present(nsamples)) then
        ns = nsamples
    else
        ns = 10000
    end if
    if (present(alpha)) then
        a = alpha
    else
        a = p05
    end if
    allocate(rst%population(ns))
    i1 = floor(half * a * ns, int32)
    i2 = ns - i1 + 1

    ! Analyze the basic data set
    rst%statistic_value = stat(x)
    rst%population(1) = rst%statistic_value

    ! Resampling Process
    call random_init(.false., .true.)
#ifdef USEOPENMP
    ! Use OpenMP to run operations in parallel
!$OMP PARALLEL DO PRIVATE(xn) SHARED(rst)
    do i = 2, ns
        ! Resample the data
        call resample(x, xn)

        ! Compute the statistic
        rst%population(i) = stat(xn)
    end do
!$OMP END PARALLEL DO
#else
    ! OpenMP is not available - run in a serial manner
    allocate(xn(n))
    do i = 2, ns
        ! Resample the data
        call resample(x, xn)

        ! Compute the statistic for the resampled data
        rst%population(i) = stat(xn)
    end do
#endif

    ! Compute the relevant quantities on the resampled statistic
    rst%upper_confidence_interval = rst%population(i2)
    rst%lower_confidence_interval = rst%population(i1)
    rst%bias = mean(rst%population) - rst%statistic_value
    rst%standard_error = standard_deviation(rst%population)
end function

! ******************************************************************************
! LINEAR REGRESSION
! ------------------------------------------------------------------------------
subroutine bootstrap_linear_least_squares(order, intercept, x, y, &
    coeffs, ymod, resid, nsamples, stats, bias, alpha, method, bscoeffs, err)
    !! Computes a linear least-squares regression to fit a set of data.
    !! Bootstrapping is utilized to gain insight into the quality of 
    !! the fit.  Resampling for the bootstrap process is a random resampling 
    !! with replacement process with the range of values limited by the 
    !! standard deviation of the original data set.
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
        !! bootstrapping analysis.  The bias is computed as the difference 
        !! between the mean of the boostrap population results for the given 
        !! parameter and the original estimate of the given parameter.
    real(real64), intent(in), optional :: alpha
        !! The significance level at which to evaluate the confidence 
        !! intervals.  The default value is 0.05 such that a 95% 
        !! confidence interval is calculated.
    procedure(bootstrap_resampling_routine), pointer, intent(in), optional :: method
        !! An optional pointer to the method to use for resampling of the data.
        !! If no method is supplied, a random resampling is utilized.
    real(real64), intent(out), optional, allocatable, target, dimension(:,:) :: bscoeffs
        !! An optional, allocatable matrix, containing the bootstrap 
        !! distributions for each parameter stored in each row of the matrix
        !! such that the resulting matrix is NCOEFFS -by- NSAMPLES.
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

    ! Local Variables
    integer(int32) :: i, j, n, ns, nc, ncoeffs, flag, nthreads, thread
    real(real64) :: alph
    real(real64), allocatable, dimension(:) :: fLocal, yLocal, rLocal
    real(real64), allocatable, target, dimension(:,:) :: coeffstorage
    real(real64), pointer, dimension(:,:) :: allcoeffs
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    procedure(bootstrap_resampling_routine), pointer :: resample
    
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
    if (present(method)) then
        resample => method
    else
        resample => random_resample
    end if
    n = size(x)
    ncoeffs = order + 1
    nc = order
    if (intercept) nc = nc + 1
    nthreads = omp_get_num_threads()

    ! Compute the fit
    call linear_least_squares(order, intercept, x, y, coeffs, &
        ymod, resid, alpha = alph, err = errmgr)
    if (errmgr%has_error_occurred()) return
    
    ! Memory Allocations
    if (present(bscoeffs)) then
        allocate(bscoeffs(ncoeffs, ns), source = zero, stat = flag)
        if (flag /= 0) then
            call report_memory_error(errmgr, "bootstrap_linear_least_squares", &
                flag)
            return
        end if
        allcoeffs => bscoeffs
    else
        allocate(coeffstorage(ncoeffs, ns), source = zero, stat = flag)
        if (flag /= 0) then
            call report_memory_error(errmgr, "bootstrap_linear_least_squares", &
                flag)
            return
        end if
        allcoeffs => coeffstorage
    end if
    allcoeffs(:,1) = coeffs

    ! Cycle over each data set and perform the fit
    call random_init(.false., .true.)
#ifdef USEOPENMP
!$OMP PARALLEL DO PRIVATE(fLocal, yLocal, rLocal) SHARED(allcoeffs)
    do i = 2, ns
        ! Get the current thread number
        ! The +1 is because OpenMP is zero-based for thread numbering
        thread = omp_get_thread_num() + 1

        ! Allocate local arrays on a per-thread basis
        if (.not.allocated(fLocal)) allocate(fLocal(n))
        if (.not.allocated(yLocal)) allocate(yLocal(n))
        if (.not.allocated(rLocal)) allocate(rLocal(n))

        ! Compute a random data set
        call resample(y, yLocal)

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
        call resample(y, yLocal)

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
            bias(i) = mean(allcoeffs(i,:)) - coeffs(i)
        end do
    end if
end subroutine

! ******************************************************************************
! NONLINEAR REGRESSION
! ------------------------------------------------------------------------------
subroutine bootstrap_nonlinear_least_squares(fun, x, y, params, ymod, resid, &
    nsamples, weights, maxp, minp, stats, alpha, controls, settings, info, &
    bias, method, bscoeffs, err)
    !! Performs a nonlinear regression to fit a model using a version
    !! of the Levenberg-Marquardt algorithm.  Bootstrapping is utilized to gain 
    !! insight into the quality of the fit.  Resampling for the bootstrap 
    !! process is a random resampling with replacement process with the 
    !! range of values limited by the standard deviation of the original 
    !! data set.
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
        !! the results of the bootstrapping analysis.  The bias is computed as 
        !! the difference between the mean of the boostrap population results 
        !! for the given parameter and the original estimate of the given 
        !! parameter.
    procedure(bootstrap_resampling_routine), pointer, intent(in), optional :: method
        !! An optional pointer to the method to use for resampling of the data.
        !! If no method is supplied, a random resampling is utilized.
    real(real64), intent(out), optional, allocatable, target, dimension(:,:) :: bscoeffs
        !! An optional, allocatable matrix, containing the bootstrap 
        !! distributions for each parameter stored in each row of the matrix
        !! such that the resulting matrix is NCOEFFS -by- NSAMPLES.
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
    
    ! Local Variables
    integer(int32) :: i, n, ns, nparams, flag
    real(real64) :: alph
    real(real64), allocatable, dimension(:) :: fLocal, yLocal, rLocal
    real(real64), allocatable, target, dimension(:,:) :: coeffstorage
    real(real64), pointer, dimension(:,:) :: allcoeffs
    procedure(bootstrap_resampling_routine), pointer :: resample
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
    if (present(method)) then
        resample => method
    else
        resample => random_resample
    end if
    n = size(x)
    nparams = size(params)

    ! Compute the fit
    call nonlinear_least_squares(fun, x, y, params, ymod, resid, &
        weights = weights, maxp = maxp, minp = minp, alpha = alph, &
        controls = controls, settings = settings, info = info, err = err)

    ! Memory Allocations
    if (present(bscoeffs)) then
        allocate(bscoeffs(nparams, ns), source = zero, stat = flag)
        if (flag /= 0) then
            call report_memory_error(errmgr, &
                "bootstrap_nonlinear_least_squares", flag)
            return
        end if
        allcoeffs => bscoeffs
    else
        allocate(coeffstorage(nparams, ns), source = zero, stat = flag)
        if (flag /= 0) then
            call report_memory_error(errmgr, &
                "bootstrap_nonlinear_least_squares", flag)
            return
        end if
        allcoeffs => coeffstorage
    end if
    allcoeffs(:,1) = params
    
    ! Define initial guesses for each step.  Base upon the results of the 
    ! initial analysis as this should provide a strong starting point for
    ! subsequent analysis
    do i = 1, nparams
        allcoeffs(i,:) = params(i)
    end do

    ! Cycle over each data set and perform the fit
    call random_init(.false., .true.)
#ifdef USEOPENMP
!$OMP PARALLEL DO
    do i = 2, ns
        ! Allocate local arrays on a per-thread basis
        if (.not.allocated(fLocal)) allocate(fLocal(n))
        if (.not.allocated(yLocal)) allocate(yLocal(n))
        if (.not.allocated(rLocal)) allocate(rLocal(n))

        ! Compute a random data set
        call resample(y, yLocal)

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
        call resample(y, yLocal)

        ! Compute the fit of the perturbed data set
        call nonlinear_least_squares(fun, x, yLocal, allcoeffs(:,i), fLocal, &
            rLocal, weights = weights, maxp = maxp, minp = minp, alpha = alph, &
            controls = controls, settings = settings, info = info)
    end do
#endif

    ! Perform the statistics calculations, if needed
    if (present(stats)) then
        ! Verify the size of stats
        if (size(stats) /= nparams) then
            call report_array_size_error(errmgr, &
                "bootstrap_nonlinear_least_squares", "stats", &
                nparams, size(stats))
            return
        end if

        ! Perform the calculations
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
            bias(i) = mean(allcoeffs(i,:)) - params(i)
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
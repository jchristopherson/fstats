module fstats_doe_models
    !! DOE model fitting and evaluation.
    use iso_fortran_env, only : int32, real64
    use fstats_errors, only : FS_INVALID_INPUT_ERROR, FS_ARRAY_SIZE_ERROR
    use fstats_regression, only : regression_statistics, covariance_matrix, &
        calculate_regression_statistics
    use fstats_experimental_design, only : doe_model
    implicit none
    private

    public :: doe_fit_model
    public :: doe_evaluate_model

    interface doe_evaluate_model
        module procedure :: doe_evaluate_model_1
        module procedure :: doe_evaluate_model_2
    end interface
contains
! ------------------------------------------------------------------------------
function doe_fit_model(nway, x, y, map, alpha) result(rst)
    use blas, only : DGEMM
    use ieee_arithmetic
    !! Fits a Taylor series model to the provided data.
    !!
    !! $$ Y = \beta_{0} + \sum_{i=1}^{n} \beta_{i} X_{i} + \sum_{i=1}^{n} 
    !! \sum_{j=1 \\ i \neq j}^{n} \beta_{ij} X_{i} X_{j} + \sum_{i=1}^{n} 
    !! \sum_{j=1}^{n} \sum_{k=1 \\ i \neq j \neq k}^{n} \beta_{ijk} X_{i} 
    !! X_{j} X_{k} + ... $$
    integer(int32), intent(in) :: nway
        !! The number of interaction levels.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! used to produce the results.
    real(real64), intent(in), dimension(:) :: y
        !! An M-element array containing the results from the M experiments.
    logical, intent(in), optional, target, dimension(:) :: map
        !! An optional array of the same size as beta that can be used to
        !! eliminate a parameter from the model (false), or keep a parameter
        !! in the model (true).  If not supplied, all parameters will be assumed
        !! to be part of the model as if the array were filled with all true
        !! values.
    real(real64), intent(in), optional :: alpha
        !! The significance level at which to evaluate the confidence 
        !! intervals.  The default value is 0.05 such that a 95% 
        !! confidence interval is calculated.
    type(doe_model) :: rst
        !! The resulting model.

    ! Local Variables
    integer(int32) :: i, j, m, n, nparam, nfactors
    logical, allocatable, target, dimension(:) :: nmap
    logical, pointer, dimension(:) :: mapptr
    real(real64) :: alph, nan
    real(real64), allocatable, dimension(:) :: coeffs, ymod, resid
    real(real64), allocatable, dimension(:,:) :: xc, c, cxt
    type(regression_statistics), allocatable, dimension(:) :: stats
    
    ! Initialization
    if (present(alpha)) then
        alph = alpha
    else
        alph = 5.0d-2
    end if
    m = size(x, 1)
    nfactors = size(x, 2)
    nan = ieee_value(nan, IEEE_QUIET_NAN)

    ! Input Checking
    if (nway < 1 .or. nway > 3) error stop FS_INVALID_INPUT_ERROR
    if (size(y) /= m) error stop FS_ARRAY_SIZE_ERROR

    ! Determine the parameter count
    nparam = 1
    if (nway >= 1) nparam = nparam + nfactors
    if (nway >= 2) nparam = nparam + nfactors * (nfactors - 1)
    if (nway >= 3) nparam = nparam + nfactors * (nfactors**2 - 1)
    
    ! Set up the map parameters
    if (present(map)) then
        if (size(map) /= nparam) error stop FS_ARRAY_SIZE_ERROR
        mapptr => map
    else
        allocate(nmap(nparam), source = .true.)
        mapptr => nmap
    end if

    ! Update the parameter count
    n = nparam
    do i = 1, nparam
        if (.not.mapptr(i)) n = n - 1
    end do
    if (n < 1) then
        error stop FS_INVALID_INPUT_ERROR
    end if

    ! Local memory allocations
    allocate(xc(m, n), c(n, n), cxt(n, m), coeffs(n))

    ! Create the design matrix
    call doe_design_matrix(nway, x, mapptr, xc)

    ! Compute the covariance matrix
    c = covariance_matrix(xc)

    ! Solve the least-squares problem (N-by-1 result)
    call DGEMM("N", "T", n, m, n, 1.0d0, c, n, xc, m, 0.0d0, cxt, n) ! C * X**T
    call DGEMM("N", "N", n, 1, m, 1.0d0, cxt, n, y, m, 0.0d0, coeffs, n) ! (C * X**T) * Y

    ! Evaluate the model and compute the residuals
    ymod = matmul(xc, coeffs)
    resid = ymod - y

    ! Estimate parameter statistics
    stats = calculate_regression_statistics(resid, coeffs, c, alph)

    ! Update output
    rst%nway = nway
    allocate(rst%coefficients(nparam))
    allocate(rst%stats(nparam))
    allocate(rst%map(nparam), source = mapptr)
    j = 0
    do i = 1, nparam
        if (mapptr(i)) then
            j = j + 1
            rst%coefficients(i) = coeffs(j)
            rst%stats(i) = stats(j)
        else
            rst%coefficients(i) = nan
            rst%stats(i)%confidence_interval = nan
            rst%stats(i)%probability = nan
            rst%stats(i)%standard_error = nan
            rst%stats(i)%t_statistic = nan
        end if
    end do
end function

! ------------------------------------------------------------------------------
subroutine doe_design_matrix(nway, x, map, c)
    !! This is an internal routine used to construct the design matrix for
    !! the DOE model of the following form:
    !!
    !! $$ Y = \beta_{0} + \sum_{i=1}^{n} \beta_{i} X_{i} + \sum_{i=1}^{n} 
    !! \sum_{j=1 \\ i \neq j}^{n} \beta_{ij} X_{i} X_{j} + \sum_{i=1}^{n} 
    !! \sum_{j=1}^{n} \sum_{k=1 \\ i \neq j \neq k}^{n} \beta_{ijk} X_{i} 
    !! X_{j} X_{k} + ... $$
    !!
    !! Up to a 3-way model is allowed.
    !!
    !! No error checking is provided.  It is assumed the arrays are sized 
    !! correctly.
    integer(int32), intent(in) :: nway
    real(real64), intent(in), dimension(:,:) :: x
    logical, intent(in), dimension(:) :: map
    real(real64), intent(out), dimension(:,:) :: c

    ! Local Variables
    integer(int32) :: i, j, k, jj, kk, m, n, np

    ! Determine the number of model parameters
    np = 0
    do i = 1, size(map)
        if (map(i)) np = np + 1
    end do

    ! Additional Initialization
    m = size(x, 1)
    n = size(x, 2)

    ! DC Term
    if (map(1)) then
        c(:,1) = 1.0d0
        jj = 2
    else
        jj = 1
    end if

    ! Main Effect
    kk = 1
    if (nway >= 1) then
        do i = 1, n
            kk = kk + 1
            if (.not.map(kk)) cycle
            c(:,jj) = x(:,i)
            jj = jj + 1
        end do
    end if

    ! Two-Way
    if (nway >= 2) then
        do i = 1, n
            do j = 1, n
                if (i == j) cycle
                kk = kk + 1
                if (.not.map(kk)) cycle
                c(:,jj) = x(:,i) * x(:,j)
                jj = jj + 1
            end do
        end do
    end if

    ! Three-Way
    if (nway >= 3) then
        do i = 1, n
            do j = 1, n
                do k = 1, n
                    if (i == j .and. j == k) cycle
                    kk = kk + 1
                    if (.not.map(kk)) cycle
                    c(:,jj) = x(:,i) * x(:,j) * x(:,k)
                    jj = jj + 1
                end do
            end do
        end do
    end if
end subroutine

! ------------------------------------------------------------------------------
function doe_evaluate_model_1(nway, beta, x, map) result(rst)
    !! Evaluates the model of the following form.
    !!
    !! $$ Y = \beta_{0} + \sum_{i=1}^{n} \beta_{i} X_{i} + \sum_{i=1}^{n} 
    !! \sum_{j=1 \\ i \neq j}^{n} \beta_{ij} X_{i} X_{j} + \sum_{i=1}^{n} 
    !! \sum_{j=1}^{n} \sum_{k=1 \\ i \neq j \neq k}^{n} \beta_{ijk} X_{i} 
    !! X_{j} X_{k} + ... $$
    integer(int32), intent(in) :: nway
        !! The number of interaction levels.  Currently, this algorithm supports
        !! a maximum of three-way interaction.
    real(real64), intent(in), dimension(:) :: beta
        !! The model coefficients.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! at which to evaluate the model.
    logical, intent(in), optional, target, dimension(:) :: map
        !! An optional array of the same size as beta that can be used to
        !! eliminate a parameter from the model (false), or keep a parameter
        !! in the model (true).  If not supplied, all parameters will be assumed
        !! to be part of the model as if the array were filled with all true
        !! values.
    real(real64), allocatable, dimension(:) :: rst
        !! The resulting M-element array.

    ! Local Variables
    integer(int32) :: m, n, nparam
    logical, pointer, dimension(:) :: mapptr
    logical, allocatable, target, dimension(:) :: nmap
    
    ! Initialization
    m = size(x, 1)
    n = size(x, 2)

    ! Input Checking
    if (nway < 1 .or. nway > 3) error stop FS_INVALID_INPUT_ERROR

    nparam = 1
    if (nway >= 1) nparam = nparam + n
    if (nway >= 2) nparam = nparam + n * (n - 1)
    if (nway >= 3) nparam = nparam + n * (n**2 - 1)
    if (size(beta) /= nparam) error stop FS_ARRAY_SIZE_ERROR

    ! Memory Allocations
    allocate(rst(m))

    ! Set up the map parameters
    if (present(map)) then
        if (size(map) /= nparam) error stop FS_ARRAY_SIZE_ERROR
        mapptr => map
    else
        allocate(nmap(nparam), source = .true.)
        mapptr => nmap
    end if

    ! Process
    call doe_eval_engine(nway, beta, x, mapptr, rst)
end function

! ------------------------------------------------------------------------------
function doe_evaluate_model_2(mdl, x) result(rst)
    !! Evaluates the model of the following form.
    !!
    !! $$ Y = \beta_{0} + \sum_{i=1}^{n} \beta_{i} X_{i} + \sum_{i=1}^{n} 
    !! \sum_{j=1 \\ i \neq j}^{n} \beta_{ij} X_{i} X_{j} + \sum_{i=1}^{n} 
    !! \sum_{j=1}^{n} \sum_{k=1 \\ i \neq j \neq k}^{n} \beta_{ijk} X_{i} 
    !! X_{j} X_{k} + ... $$
    class(doe_model), intent(in) :: mdl
        !! The model to evaluate.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! at which to evaluate the model.
    real(real64), allocatable, dimension(:) :: rst
        !! The resulting M-element array.

    ! Process
    rst = doe_evaluate_model_1(mdl%nway, mdl%coefficients, x, mdl%map)
end function

! ------------------------------------------------------------------------------
subroutine doe_eval_engine(nway, beta, x, map, y)
    ! Driver routine for "doe_evaluate_model" that performs the actual 
    ! calculations but forgoes any error checking.  This should not be exposed 
    ! as part of the public API.
    integer(int32), intent(in) :: nway
    real(real64), intent(in), dimension(:) :: beta
    real(real64), intent(in), dimension(:,:) :: x
    logical, intent(in), dimension(:) :: map
    real(real64), intent(out), dimension(:) :: y

    ! Local Variables
    integer(int32) :: i1, i2, n

    ! Initialization
    n = size(x, 2)
    if (map(1)) then
        y = beta(1)
    else
        y = 0.0d0
    end if

    ! Process
    if (nway >= 1) then
        i1 = 2
        i2 = i1 + n - 1
        call doe_eval_1(beta(i1:i2), x, map(i1:i2), y)
    end if
    if (nway >= 2) then
        i1 = i2 + 1
        i2 = i1 + n * (n - 1) - 1
        call doe_eval_2(beta(i1:i2), x, map(i1:i2), y)
    end if
    if (nway >= 3) then
        i1 = i2 + 1
        i2 = i1 + n * (n**2 - 1) - 1
        call doe_eval_3(beta(i1:i2), x, map(i1:i2), y)
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine doe_eval_1(beta, x, map, y)
    !! Evaluates the main effect term.
    !!
    !! $$ Y = Y + /sum_{i=1}^{n} \beta_{i} X_{i} $$
    real(real64), intent(in), dimension(:) :: beta
        !! The model coefficients for just this portion of the model.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! at which to evaluate the model.
    logical, intent(in), dimension(:) :: map
        !! The usage map corresponding to the model coefficients for just this
        !! portion of the model.
    real(real64), intent(inout), dimension(:) :: y
        !! On input, an M-element array containing the existing portion of the 
        !! model.  On output, this array is updated to include the main effects.

    ! Local Variables
    integer(int32) :: i, n

    ! Initialization
    n = size(x, 2)

    ! Process
    do i = 1, n
        if (.not.map(i)) cycle
        y = y + beta(i) * x(:,i)
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine doe_eval_2(beta, x, map, y)
    !! Evaluates the two-way interaction term.
    !!
    !! $$ Y = Y + /sum_{i=1}^{n} /sum_{j=1 // i /neq j}^{n} \beta_{i} X_{i} 
    !! X_{j} $$
    real(real64), intent(in), dimension(:) :: beta
        !! The model coefficients for just this portion of the model.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! at which to evaluate the model.
    logical, intent(in), dimension(:) :: map
        !! The usage map corresponding to the model coefficients for just this
        !! portion of the model.
    real(real64), intent(inout), dimension(:) :: y
        !! On input, an M-element array containing the existing portion of the 
        !! model.  On output, this array is updated to include the two-way
        !! interactions.

    ! Local Variables
    integer(int32) :: i, j, k, n

    ! Initialization
    n = size(x, 2)

    ! Process
    k = 0
    do i = 1, n
        do j = 1, n
            if (i == j) cycle
            k = k + 1
            if (.not.map(k)) cycle
            y = y + beta(k) * x(:,i) * x(:,j)
        end do
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine doe_eval_3(beta, x, map, y)
    !! Evaluates the three-way interaction term.
    !!
    !! $$ Y = Y + /sum_{i=1}^{n} /sum_{j=1 // i /neq j}^{n} /sum_{k=1 // i /neq j /neq k}^{n} \beta_{i} X_{i} 
    !! X_{j} X_{k} $$
    real(real64), intent(in), dimension(:) :: beta
        !! The model coefficients for just this portion of the model.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! at which to evaluate the model.
    logical, intent(in), dimension(:) :: map
        !! The usage map corresponding to the model coefficients for just this
        !! portion of the model.
    real(real64), intent(inout), dimension(:) :: y
        !! On input, an M-element array containing the existing portion of the 
        !! model.  On output, this array is updated to include the three-way
        !! interactions.

    ! Local Variables
    integer(int32) :: i, j, k, ii, n

    ! Initialization
    n = size(x, 2)

    ! Process
    ii = 0
    do i = 1, n
        do j = 1, n
            do k = 1, n
                if (i == j .and. j == k) cycle
                ii = ii + 1
                if (.not.map(ii)) cycle
                y = y + beta(ii) * x(:,i) * x(:,j) * x(:,k)
            end do
        end do
    end do
end subroutine

end module fstats_doe_models

submodule (fstats) regression_implementation
    use linalg
    use fstats_errors
    use external
    implicit none
contains
! ------------------------------------------------------------------------------
module subroutine coefficient_matrix_real64(order, intercept, x, c, err)
    ! Arguments
    integer(int32), intent(in) :: order
    logical, intent(in) :: intercept
    real(real64), intent(in) :: x(:)
    real(real64), intent(out) :: c(:,:)
    class(errors), intent(inout), optional, target :: err

    ! Parameters
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, start, npts, ncols
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    npts = size(x)
    ncols = order
    if (intercept) ncols = ncols + 1

    ! Input Check
    if (order < 1) then
        call errmgr%report_error("coefficient_matrix_real64", &
            "The model order must be at least one.", FS_INVALID_INPUT_ERROR)
        return
    end if
    if (size(c, 1) /= npts .or. size(c, 2) /= ncols) then
        call report_matrix_size_error(errmgr, "coefficient_matrix_real64", &
            "c", npts, ncols, size(c, 1), size(c, 2))
        return
    end if

    ! Process
    if (intercept) then
        c(:,1) = one
        c(:,2) = x
        start = 3
    else
        c(:,1) = x
        start = 2
    end if
    if (start >= ncols) return
    do i = start, ncols
        c(:,i) = c(:,i-1) * x
    end do
end subroutine

! --------------------
module subroutine coefficient_matrix_real32(order, intercept, x, c, err)
    ! Arguments
    integer(int32), intent(in) :: order
    logical, intent(in) :: intercept
    real(real32), intent(in) :: x(:)
    real(real32), intent(out) :: c(:,:)
    class(errors), intent(inout), optional, target :: err

    ! Parameters
    real(real32), parameter :: one = 1.0

    ! Local Variables
    integer(int32) :: i, start, npts, ncols
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    npts = size(x)
    ncols = order
    if (intercept) ncols = ncols + 1

    ! Input Check
    if (order < 1) then
        call errmgr%report_error("coefficient_matrix_real32", &
            "The model order must be at least one.", FS_INVALID_INPUT_ERROR)
        return
    end if
    if (size(c, 1) /= npts .or. size(c, 2) /= ncols) then
        call report_matrix_size_error(errmgr, "coefficient_matrix_real32", &
            "c", npts, ncols, size(c, 1), size(c, 2))
        return
    end if

    ! Process
    if (intercept) then
        c(:,1) = one
        c(:,2) = x
        start = 3
    else
        c(:,1) = x
        start = 2
    end if
    if (start >= ncols) return
    do i = start, ncols
        c(:,i) = c(:,i-1) * x
    end do
end subroutine

! ------------------------------------------------------------------------------
module subroutine covariance_matrix_real64(x, c, err)
    ! Arguments
    real(real64), intent(in) :: x(:,:)
    real(real64), intent(out) :: c(:,:)
    class(errors), intent(inout), optional, target :: err

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    integer(int32) :: npts, ncoeffs, flag
    real(real64), allocatable :: xtx(:,:)
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    npts = size(x, 1)
    ncoeffs = size(x, 2)

    ! Input Checking
    if (size(c, 1) /= ncoeffs .or. size(c, 2) /= ncoeffs) then
        call report_matrix_size_error(errmgr, "covariance_matrix_real64", &
            "c", ncoeffs, ncoeffs, size(c, 1), size(c, 2))
        return
    end if

    ! Local Memory Allocation
    allocate(xtx(ncoeffs, ncoeffs), stat = flag)
    if (flag /= 0) then
        call report_memory_error(errmgr, "covariance_matrix_real64", flag)
        return
    end if

    ! Compute X**T * X
    call DGEMM("T", "N", ncoeffs, ncoeffs, npts, one, x, npts, x, npts, &
        zero, xtx, ncoeffs)
    
    ! Compute the inverse of X**T * X to obtain the covariance matrix
    call mtx_pinverse(xtx, c, err = errmgr)
    if (errmgr%has_error_occurred()) return
end subroutine

! --------------------
module subroutine covariance_matrix_real32(x, c, err)
    ! Arguments
    real(real32), intent(in) :: x(:,:)
    real(real32), intent(out) :: c(:,:)
    class(errors), intent(inout), optional, target :: err

    ! Parameters
    real(real32), parameter :: zero = 0.0d0
    real(real32), parameter :: one = 1.0d0

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    integer(int32) :: npts, ncoeffs, flag
    real(real32), allocatable :: xtx(:,:)
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    npts = size(x, 1)
    ncoeffs = size(x, 2)

    ! Input Checking
    if (size(c, 1) /= ncoeffs .or. size(c, 2) /= ncoeffs) then
        call report_matrix_size_error(errmgr, "covariance_matrix_real32", &
            "c", ncoeffs, ncoeffs, size(c, 1), size(c, 2))
        return
    end if

    ! Local Memory Allocation
    allocate(xtx(ncoeffs, ncoeffs), stat = flag)
    if (flag /= 0) then
        call report_memory_error(errmgr, "covariance_matrix_real64", flag)
        return
    end if

    ! Compute X**T * X
    call SGEMM("T", "N", ncoeffs, ncoeffs, npts, one, x, npts, x, npts, &
        zero, xtx, ncoeffs)
    
    ! Compute the inverse of X**T * X to obtain the covariance matrix
    call r32_inverse(xtx, c, errmgr)
    if (errmgr%has_error_occurred()) return
end subroutine

! ------------------------------------------------------------------------------
module subroutine linear_least_squares_real64(order, intercept, x, y, coeffs, &
    ymod, resid, stats, alpha, err)
    ! Arguments
    integer(int32), intent(in) :: order
    logical, intent(in) :: intercept
    real(real64), intent(in) :: x(:), y(:)
    real(real64), intent(out) :: coeffs(:), ymod(:), resid(:)
    type(regression_statistics), intent(out), optional :: stats(:)
    real(real64), intent(in), optional :: alpha
    class(errors), intent(inout), optional, target :: err

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: half = 0.5d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, npts, ncols, ncoeffs, flag
    real(real64) :: alph, var, df, ssr, talpha
    real(real64), allocatable :: a(:,:), c(:,:), cxt(:,:)
    type(t_distribution) :: dist
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    npts = size(x)
    ncoeffs = order + 1
    ncols = order
    if (intercept) ncols = ncols + 1
    alph = 0.05d0
    if (present(alpha)) alph = alpha

    ! Input Check
    if (order < 1) then
        call errmgr%report_error("linear_least_squares_real64", &
            "The model order must be at least one.", FS_INVALID_INPUT_ERROR)
        return
    end if
    if (size(y) /= npts) then
        call report_array_size_error(errmgr, "linear_least_squares_real64", &
            "y", npts, size(y))
        return
    end if
    if (size(coeffs) /= ncoeffs) then
        call report_array_size_error(errmgr, "linear_least_squares_real64", &
            "coeffs", ncoeffs, size(coeffs))
        return
    end if
    if (size(ymod) /= npts) then
        call report_array_size_error(errmgr, "linear_least_squares_real64", &
            "ymod", npts, size(ymod))
        return
    end if
    if (size(resid) /= npts) then
        call report_array_size_error(errmgr, "linear_least_squares_real64", &
            "resid", npts, size(resid))
        return
    end if
    if (present(stats)) then
        if (size(stats) /= ncols) then
            call report_array_size_error(errmgr, &
                "linear_least_squares_real64", "stats", ncols, size(stats))
            return
        end if
    end if

    ! Memory Allocation
    allocate(a(npts, ncols), stat = flag)
    if (flag == 0) allocate(c(ncols, ncols), stat = flag)
    if (flag == 0) allocate(cxt(ncols, npts), stat = flag)
    if (flag /= 0) then
        call report_memory_error(errmgr, "linear_least_squares_real64", flag)
        return
    end if

    ! Compute the coefficient matrix
    call coefficient_matrix(order, intercept, x, a, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute the covariance matrix
    call covariance_matrix(a, c, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute the coefficients (NCOLS-by-1)
    call DGEMM("N", "T", ncols, npts, ncols, one, c, ncols, a, npts, zero, &
        cxt, ncols)     ! C * X**T

    i = 2
    coeffs(1) = zero
    if (intercept) i = 1
    call DGEMM("N", "N", ncols, 1, npts, one, cxt, ncols, y, npts, zero, &
        coeffs(i:ncoeffs), ncols)  ! (C * X**T) * Y

    ! Evaluate the model and compute the residuals
    call DGEMM("N", "N", npts, 1, ncols, one, a, npts, coeffs(i:ncoeffs), &
        ncols, zero, ymod, npts)
    resid = ymod - y

    ! If the user doesn't want the statistics calculations we can stop now
    if (.not.present(stats)) return
    
    ! Start the process of computing statistics
    stats = calculate_regression_statistics(resid, coeffs(1:ncols), c, alph, &
        errmgr)
end subroutine

! --------------------
module subroutine linear_least_squares_real32(order, intercept, x, y, coeffs, &
    ymod, resid, stats, alpha, err)
    ! Arguments
    integer(int32), intent(in) :: order
    logical, intent(in) :: intercept
    real(real32), intent(in) :: x(:), y(:)
    real(real32), intent(out) :: coeffs(:), ymod(:), resid(:)
    type(regression_statistics), intent(out), optional :: stats(:)
    real(real32), intent(in), optional :: alpha
    class(errors), intent(inout), optional, target :: err

    ! Parameters
    real(real32), parameter :: zero = 0.0
    real(real32), parameter :: half = 0.5
    real(real32), parameter :: one = 1.0

    ! Local Variables
    integer(int32) :: i, npts, ncols, ncoeffs, flag
    real(real32) :: alph, var, df, ssr, talpha
    real(real32), allocatable :: a(:,:), c(:,:), cxt(:,:)
    type(t_distribution) :: dist
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    npts = size(x)
    ncoeffs = order + 1
    ncols = order
    if (intercept) ncols = ncols + 1
    alph = 0.05
    if (present(alpha)) alph = alpha

    ! Input Check
    if (order < 1) then
        call errmgr%report_error("linear_least_squares_real32", &
            "The model order must be at least one.", FS_INVALID_INPUT_ERROR)
        return
    end if
    if (size(y) /= npts) then
        call report_array_size_error(errmgr, "linear_least_squares_real32", &
            "y", npts, size(y))
        return
    end if
    if (size(coeffs) /= ncoeffs) then
        call report_array_size_error(errmgr, "linear_least_squares_real32", &
            "coeffs", ncoeffs, size(coeffs))
        return
    end if
    if (size(ymod) /= npts) then
        call report_array_size_error(errmgr, "linear_least_squares_real32", &
            "ymod", npts, size(ymod))
        return
    end if
    if (size(resid) /= npts) then
        call report_array_size_error(errmgr, "linear_least_squares_real32", &
            "resid", npts, size(resid))
        return
    end if
    if (present(stats)) then
        if (size(stats) /= ncols) then
            call report_array_size_error(errmgr, &
                "linear_least_squares_real32", "stats", ncols, size(stats))
            return
        end if
    end if

    ! Memory Allocation
    allocate(a(npts, ncols), stat = flag)
    if (flag == 0) allocate(c(ncols, ncols), stat = flag)
    if (flag == 0) allocate(cxt(ncols, npts), stat = flag)
    if (flag /= 0) then
        call report_memory_error(errmgr, "linear_least_squares_real32", flag)
        return
    end if

    ! Compute the coefficient matrix
    call coefficient_matrix(order, intercept, x, a, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute the covariance matrix
    call covariance_matrix(a, c, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute the coefficients (NCOLS-by-1)
    call SGEMM("N", "T", ncols, npts, ncols, one, c, ncols, a, npts, zero, &
        cxt, ncols)     ! C * X**T

    i = 2
    coeffs(1) = zero
    if (intercept) i = 1
    call SGEMM("N", "N", ncols, 1, npts, one, cxt, ncols, y, npts, zero, &
        coeffs(i:ncoeffs), ncols)  ! (C * X**T) * Y

    ! Evaluate the model and compute the residuals
    call SGEMM("N", "N", npts, 1, ncols, one, a, npts, coeffs(i:ncoeffs), &
        ncols, zero, ymod, npts)
    resid = ymod - y

    ! If the user doesn't want the statistics calculations we can stop now
    if (.not.present(stats)) return
    
    ! Start the process of computing statistics
    stats = calculate_regression_statistics(resid, coeffs(1:ncols), c, alph, &
        errmgr)
end subroutine

! ------------------------------------------------------------------------------
module function calculate_regression_stats_r64(resid, params, c, alpha, err) &
    result(rst)
    ! Arguments
    real(real64), intent(in) :: resid(:), params(:), c(:,:)
    real(real64), intent(in), optional :: alpha
    class(errors), intent(inout), optional, target :: err
    type(regression_statistics), allocatable :: rst(:)

    ! Parameters
    real(real64), parameter :: p05 = 0.05d0
    real(real64), parameter :: half = 0.5d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, m, n, dof, flag
    real(real64) :: a, ssr, var, talpha
    type(t_distribution) :: dist
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Initialization
    m = size(resid)
    n = size(params)
    dof = m - n
    if (present(alpha)) then
        a = alpha
    else
        a = p05
    end if
    allocate(rst(n), stat = flag)
    if (flag /= 0) then
    end if

    ! Input Checking
    if (size(c, 1) /= n .or. size(c, 2) /= n) then
    end if

    ! Process
    ssr = norm2(resid)**2   ! sum of the squares of the residual
    var = ssr / dof
    dist%dof = real(dof, real64)
    talpha = confidence_interval(dist, a, one, 1)
    do i = 1, n
        rst(i)%standard_error = sqrt(var * c(i,i))
        rst(i)%t_statistic = params(i) / rst(i)%standard_error
        rst(i)%probability = regularized_beta( &
            half * dof, &
            half, &
            real(dof, real64) / (dof + (rst(i)%t_statistic)**2) &
        )
        rst(i)%confidence_interval = talpha * rst(i)%standard_error
    end do
end function

! --------------------
module function calculate_regression_stats_r32(resid, params, c, alpha, err) &
    result(rst)
    ! Arguments
    real(real32), intent(in) :: resid(:), params(:), c(:,:)
    real(real32), intent(in), optional :: alpha
    class(errors), intent(inout), optional, target :: err
    type(regression_statistics), allocatable :: rst(:)

    ! Parameters
    real(real32), parameter :: p05 = 0.05
    real(real64), parameter :: half = 0.5d0
    real(real32), parameter :: one = 1.0

    ! Local Variables
    integer(int32) :: i, m, n, dof, flag
    real(real32) :: a, ssr, var, talpha
    type(t_distribution) :: dist
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Initialization
    m = size(resid)
    n = size(params)
    dof = m - n
    if (present(alpha)) then
        a = alpha
    else
        a = p05
    end if
    allocate(rst(n), stat = flag)
    if (flag /= 0) then
    end if

    ! Input Checking
    if (size(c, 1) /= n .or. size(c, 2) /= n) then
    end if

    ! Process
    ssr = norm2(resid)**2   ! sum of the squares of the residual
    var = ssr / dof
    dist%dof = real(dof, real32)
    talpha = confidence_interval(dist, a, one, 1)
    do i = 1, n
        rst(i)%standard_error = sqrt(var * c(i,i))
        rst(i)%t_statistic = params(i) / rst(i)%standard_error
        rst(i)%probability = regularized_beta( &
            half * dof, &
            half, &
            real(dof, real64) / (dof + (rst(i)%t_statistic)**2) &
        )
        rst(i)%confidence_interval = talpha * rst(i)%standard_error
    end do
end function

! ******************************************************************************
! PRIVATE ROUTINES
! ------------------------------------------------------------------------------
subroutine r32_inverse(x, xinv, err)
    ! Arguments
    real(real32), intent(in) :: x(:,:)
    real(real32), intent(out) :: xinv(:,:)
    class(errors), intent(inout) :: err

    ! Variables
    integer(int32) :: m, n
    real(real64), allocatable :: xc(:,:), xcinv(:,:)

    ! Process
    m = size(x, 1)
    n = size(x, 2)
    allocate(xcinv(n, m))
    xc = real(x, real64)
    call mtx_pinverse(xc, xcinv, err = err)
    xinv = real(xcinv, real32)
end subroutine

! ------------------------------------------------------------------------------
end submodule
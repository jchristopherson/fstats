submodule (fstats) regression_implementation
    use linalg
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
    character(len = 256) :: errmsg
    
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
        write(errmsg, 100) "The output matrix was expected to by ", npts, &
            "-by-", ncols, ", but was found to be ", size(c, 1), "-by-", &
            size(c, 2), "."
        call errmgr%report_error("coefficient_matrix_real64", trim(errmsg), &
            FS_ARRAY_SIZE_ERROR)
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

    ! Formatting
100 format(A, I0, A, I0, A, I0, A, I0, A)
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
    character(len = 256) :: errmsg
    
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
        write(errmsg, 100) "The output matrix was expected to by ", npts, &
            "-by-", ncols, ", but was found to be ", size(c, 1), "-by-", &
            size(c, 2), "."
        call errmgr%report_error("coefficient_matrix_real32", trim(errmsg), &
            FS_ARRAY_SIZE_ERROR)
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

    ! Formatting
100 format(A, I0, A, I0, A, I0, A, I0, A)
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
    character(len = 256) :: errmsg
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
        write(errmsg, 100) "The covariance matrix was expected to be ", &
            ncoeffs, "-by-", ncoeffs, ", but was found to be ", size(c, 1), &
            "-by-", size(c, 2), "."
        call errmgr%report_error("covariance_matrix_real64", trim(errmsg), &
            FS_ARRAY_SIZE_ERROR)
        return
    end if

    ! Local Memory Allocation
    allocate(xtx(ncoeffs, ncoeffs), stat = flag)
    if (flag /= 0) then
        write(errmsg, 101) "Memory allocation error code ", flag, "."
        call errmgr%report_error("covariance_matrix_real64", &
            trim(errmsg), FS_OUT_OF_MEMORY_ERROR)
        return
    end if

    ! Compute X**T * X
    call DGEMM("T", "N", ncoeffs, ncoeffs, npts, one, x, npts, x, npts, &
        zero, xtx, ncoeffs)
    
    ! Compute the inverse of X**T * X to obtain the covariance matrix
    call mtx_pinverse(xtx, c, err = errmgr)
    if (errmgr%has_error_occurred()) return

    ! Formatting
100 format(A, I0, A, I0, A, I0, A, I0, A)
101 format(A, I0, A)
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
    character(len = 256) :: errmsg
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
        write(errmsg, 100) "The covariance matrix was expected to be ", &
            ncoeffs, "-by-", ncoeffs, ", but was found to be ", size(c, 1), &
            "-by-", size(c, 2), "."
        call errmgr%report_error("covariance_matrix_real32", trim(errmsg), &
            FS_ARRAY_SIZE_ERROR)
        return
    end if

    ! Local Memory Allocation
    allocate(xtx(ncoeffs, ncoeffs), stat = flag)
    if (flag /= 0) then
        write(errmsg, 101) "Memory allocation error code ", flag, "."
        call errmgr%report_error("covariance_matrix_real32", &
            trim(errmsg), FS_OUT_OF_MEMORY_ERROR)
        return
    end if

    ! Compute X**T * X
    call SGEMM("T", "N", ncoeffs, ncoeffs, npts, one, x, npts, x, npts, &
        zero, xtx, ncoeffs)
    
    ! Compute the inverse of X**T * X to obtain the covariance matrix
    call r32_inverse(xtx, c, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Formatting
100 format(A, I0, A, I0, A, I0, A, I0, A)
101 format(A, I0, A)
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
    character(len = 256) :: errmsg
    
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
        write(errmsg, 100) &
            "The dependent data array was expected to be of length ", npts, &
            ", but was found to be of length ", size(y), "."
        call errmgr%report_error("linear_least_squares_real64", &
            trim(errmsg), FS_ARRAY_SIZE_ERROR)
    end if
    if (size(coeffs) /= ncoeffs) then
        write(errmsg, 100) &
            "The coefficient array was expected to be of length ", ncoeffs, &
            ", but was found to be of length ", size(coeffs), "."
        call errmgr%report_error("linear_least_squares_real64", &
            trim(errmsg), FS_ARRAY_SIZE_ERROR)
    end if
    if (size(ymod) /= npts) then
        write(errmsg, 100) &
            "The modeled data array was expected to be of length ", npts, &
            ", but was found to be of length ", size(ymod), "."
        call errmgr%report_error("linear_least_squares_real64", &
            trim(errmsg), FS_ARRAY_SIZE_ERROR)
    end if
    if (size(resid) /= npts) then
        write(errmsg, 100) &
            "The residuals array was expected to be of length ", npts, &
            ", but was found to be of length ", size(resid), "."
        call errmgr%report_error("linear_least_squares_real64", &
            trim(errmsg), FS_ARRAY_SIZE_ERROR)
    end if
    if (present(stats)) then
        if (size(stats) /= ncols) then
            write(errmsg, 100) &
                "The statistics array was expected to be of length ", ncols, &
                ", but was found to be of length ", size(stats), "."
            call errmgr%report_error("linear_least_squares_real64", &
                trim(errmsg), FS_ARRAY_SIZE_ERROR)
        end if
    end if

    ! Memory Allocation
    allocate(a(npts, ncols), stat = flag)
    if (flag == 0) allocate(c(ncols, ncols), stat = flag)
    if (flag == 0) allocate(cxt(ncols, npts), stat = flag)
    if (flag /= 0) then
        write(errmsg, 101) "Memory allocation error code ", flag, "."
        call errmgr%report_error("linear_least_squares_real64", &
            trim(errmsg), FS_OUT_OF_MEMORY_ERROR)
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
    ! df = npts - ncols
    ! ssr = norm2(resid)**2   ! sum of the squares of the residual
    ! var = ssr / df
    ! dist%dof = df
    ! talpha = confidence_interval(dist, alph, one, 1)
    ! do i = 1, ncols
    !     stats(i)%standard_error = sqrt(var * c(i,i))
    !     stats(i)%t_statistic = coeffs(i) / stats(i)%standard_error
    !     stats(i)%probability = regularized_beta( &
    !         half * df, &
    !         half, &
    !         df / (df + (stats(i)%t_statistic)**2) &
    !     )
    !     stats(i)%confidence_interval = talpha * stats(i)%standard_error
    ! end do
    stats = calculate_regression_statistics(resid, coeffs(1:ncols), c, alph, &
        errmgr)

    ! Formatting
100 format(A, I0, A, I0, A)
101 format(A, I0, A)
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
    character(len = 256) :: errmsg
    
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
        write(errmsg, 100) &
            "The dependent data array was expected to be of length ", npts, &
            ", but was found to be of length ", size(y), "."
        call errmgr%report_error("linear_least_squares_real32", &
            trim(errmsg), FS_ARRAY_SIZE_ERROR)
    end if
    if (size(coeffs) /= ncoeffs) then
        write(errmsg, 100) &
            "The coefficient array was expected to be of length ", ncoeffs, &
            ", but was found to be of length ", size(coeffs), "."
        call errmgr%report_error("linear_least_squares_real32", &
            trim(errmsg), FS_ARRAY_SIZE_ERROR)
    end if
    if (size(ymod) /= npts) then
        write(errmsg, 100) &
            "The modeled data array was expected to be of length ", npts, &
            ", but was found to be of length ", size(ymod), "."
        call errmgr%report_error("linear_least_squares_real32", &
            trim(errmsg), FS_ARRAY_SIZE_ERROR)
    end if
    if (size(resid) /= npts) then
        write(errmsg, 100) &
            "The residuals array was expected to be of length ", npts, &
            ", but was found to be of length ", size(resid), "."
        call errmgr%report_error("linear_least_squares_real32", &
            trim(errmsg), FS_ARRAY_SIZE_ERROR)
    end if
    if (present(stats)) then
        if (size(stats) /= ncols) then
            write(errmsg, 100) &
                "The statistics array was expected to be of length ", ncols, &
                ", but was found to be of length ", size(stats), "."
            call errmgr%report_error("linear_least_squares_real32", &
                trim(errmsg), FS_ARRAY_SIZE_ERROR)
        end if
    end if

    ! Memory Allocation
    allocate(a(npts, ncols), stat = flag)
    if (flag == 0) allocate(c(ncols, ncols), stat = flag)
    if (flag == 0) allocate(cxt(ncols, npts), stat = flag)
    if (flag /= 0) then
        write(errmsg, 101) "Memory allocation error code ", flag, "."
        call errmgr%report_error("linear_least_squares_real32", &
            trim(errmsg), FS_OUT_OF_MEMORY_ERROR)
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
    ! df = npts - ncols
    ! ssr = norm2(resid)**2   ! sum of the squares of the residual
    ! var = ssr / df
    ! dist%dof = df
    ! talpha = confidence_interval(dist, alph, one, 1)
    ! do i = 1, ncols
    !     stats(i)%standard_error = sqrt(var * c(i,i))
    !     stats(i)%t_statistic = coeffs(i) / stats(i)%standard_error
    !     stats(i)%probability = regularized_beta( &
    !         half * df, &
    !         half, &
    !         real(df / (df + (stats(i)%t_statistic)**2), real32) &
    !     )
    !     stats(i)%confidence_interval = talpha * stats(i)%standard_error
    ! end do
    stats = calculate_regression_statistics(resid, coeffs(1:ncols), c, alph, &
        errmgr)

    ! Formatting
100 format(A, I0, A, I0, A)
101 format(A, I0, A)
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
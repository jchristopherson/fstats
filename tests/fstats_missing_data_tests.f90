module fstats_missing_data_tests
    use iso_fortran_env
    use fortran_test_helper
    use fstats
    implicit none

contains
! ------------------------------------------------------------------------------
function test_missing_value() result(rst)
    logical :: rst

    ! Local Variables
    real(real64) :: x(5)
    logical :: flags(5)

    ! Initialization
    rst = .true.
    x = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0]
    x(3) = missing_value()

    ! Test 1 - the flag must register as missing
    if (.not.is_missing(missing_value())) then
        rst = .false.
        print "(A)", "TEST FAILED: test_missing_value -1"
    end if

    ! Test 2 - ordinary values must not register as missing
    if (is_missing(0.0d0) .or. is_missing(-1.0d3)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_missing_value -2"
    end if

    ! Test 3 - the routine must operate elementally
    flags = is_missing(x)
    if (count(flags) /= 1 .or. .not.flags(3)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_missing_value -3"
    end if
end function

! ------------------------------------------------------------------------------
function test_knn_impute() result(rst)
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-10

    ! Local Variables
    integer(int32) :: i
    real(real64) :: x(6,2), xc(6,2), xw(6,2), ans

    ! Initialization - two well-separated clusters with a single missing value
    ! in the second cluster
    rst = .true.
    x(:,1) = [0.0d0, 1.0d-1, 2.0d-1, 5.0d0, 5.1d0, 5.2d0]
    x(:,2) = [1.0d1, 1.01d1, 1.02d1, 5.0d1, 5.01d1, 0.0d0]
    x(6,2) = missing_value()

    call knn_impute(x, xc, k = 2, weighted = .false.)

    ! Test 1 - no missing values may remain
    if (any(is_missing(xc))) then
        rst = .false.
        print "(A)", "TEST FAILED: test_knn_impute -1"
    end if

    ! Test 2 - the observed values must be left untouched
    do i = 1, 5
        if (.not.assert(x(i,:), xc(i,:), tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_knn_impute -2"
            exit
        end if
    end do
    if (.not.assert(x(6,1), xc(6,1), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_knn_impute -3"
    end if

    ! Test 3 - the two nearest donors are rows 5 and 4; the unweighted estimate
    ! is therefore the mean of their second variable
    ans = 0.5d0 * (x(5,2) + x(4,2))
    if (.not.assert(ans, xc(6,2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_knn_impute -4"
    end if

    ! Test 4 - the distance-weighted estimate must lie between the two donor
    ! values, and be nearer the closer of the two
    call knn_impute(x, xw, k = 2, weighted = .true.)
    if (xw(6,2) <= x(4,2) .or. xw(6,2) >= x(5,2) .or. xw(6,2) <= ans) then
        rst = .false.
        print "(A)", "TEST FAILED: test_knn_impute -5"
    end if
end function

! ------------------------------------------------------------------------------
function test_em_impute_complete() result(rst)
    !! With no missing data the EM estimates must reduce to the usual
    !! maximum-likelihood estimates, and the data must be returned unaltered.
    logical :: rst

    ! Parameters
    integer(int32), parameter :: n = 50
    integer(int32), parameter :: m = 3
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    integer(int32) :: i, j
    real(real64) :: x(n,m), xc(n,m), mu(m), sigma(m,m), ans(m,m), amu(m)

    ! Initialization
    rst = .true.
    call random_number(x)
    call em_impute(x, xc, mu = mu, sigma = sigma)

    ! Test 1 - the data must be returned unaltered
    if (.not.assert(x, xc, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_em_impute_complete -1"
    end if

    ! Test 2 - the mean vector
    do j = 1, m
        amu(j) = sum(x(:,j)) / real(n, real64)
    end do
    if (.not.assert(amu, mu, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_em_impute_complete -2"
    end if

    ! Test 3 - the covariance matrix, noting that the maximum-likelihood
    ! estimator normalizes by N rather than N - 1
    ans = 0.0d0
    do i = 1, n
        do j = 1, m
            ans(:,j) = ans(:,j) + (x(i,:) - amu) * (x(i,j) - amu(j))
        end do
    end do
    ans = ans / real(n, real64)
    if (.not.assert(ans, sigma, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_em_impute_complete -3"
    end if
end function

! ------------------------------------------------------------------------------
function test_em_impute() result(rst)
    logical :: rst

    ! Parameters
    integer(int32), parameter :: n = 2000
    integer(int32), parameter :: m = 2
    real(real64), parameter :: tol = 1.5d-1

    ! Local Variables
    integer(int32) :: i, j, niter
    real(real64) :: x(n,m), xc(n,m), mu(m), sigma(m,m), u
    real(real64) :: mtrue(m), l(m,m)

    ! Initialization - draw a correlated, bivariate normal sample
    rst = .true.
    mtrue = [1.0d0, -2.0d0]
    l = reshape([1.0d0, 8.0d-1, 0.0d0, 6.0d-1], [m, m])
    do i = 1, n
        x(i,:) = mtrue + matmul(l, box_muller_sample(0.0d0, 1.0d0, m))
    end do

    ! Remove roughly 20% of the values from the second variable only such that
    ! the mechanism is missing at random
    do i = 1, n
        call random_number(u)
        if (u < 2.0d-1) x(i,2) = missing_value()
    end do

    call em_impute(x, xc, mu = mu, sigma = sigma, niter = niter)

    ! Test 1 - no missing values may remain
    if (any(is_missing(xc))) then
        rst = .false.
        print "(A)", "TEST FAILED: test_em_impute -1"
    end if

    ! Test 2 - the observed values must be left untouched
    do i = 1, n
        do j = 1, m
            if (is_missing(x(i,j))) cycle
            if (.not.assert(x(i,j), xc(i,j), 1.0d-10)) then
                rst = .false.
                print "(A)", "TEST FAILED: test_em_impute -2"
                exit
            end if
        end do
    end do

    ! Test 3 - the mean vector must be recovered
    if (.not.assert(mtrue, mu, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_em_impute -3"
    end if

    ! Test 4 - the covariance matrix must be recovered
    if (.not.assert(matmul(l, transpose(l)), sigma, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_em_impute -4"
    end if

    ! Test 5 - the iteration must have converged
    if (niter < 1 .or. niter >= 500) then
        rst = .false.
        print "(A)", "TEST FAILED: test_em_impute -5"
    end if
end function

! ------------------------------------------------------------------------------
function test_multiple_impute() result(rst)
    logical :: rst

    ! Parameters
    integer(int32), parameter :: n = 500
    integer(int32), parameter :: m = 2
    integer(int32), parameter :: nk = 5
    real(real64), parameter :: tol = 2.5d-1

    ! Parameters
    integer(int32) :: i, j, k
    logical :: varies
    real(real64) :: x(n,m), xc(n,m,nk), mtrue(m), l(m,m), u, avg(m)

    ! Initialization
    rst = .true.
    mtrue = [1.0d0, -2.0d0]
    l = reshape([1.0d0, 8.0d-1, 0.0d0, 6.0d-1], [m, m])
    do i = 1, n
        x(i,:) = mtrue + matmul(l, box_muller_sample(0.0d0, 1.0d0, m))
    end do
    do i = 1, n
        call random_number(u)
        if (u < 2.0d-1) x(i,2) = missing_value()
    end do

    call multiple_impute(x, xc, nburn = 10)

    ! Test 1 - no missing values may remain in any of the imputations
    if (any(is_missing(xc))) then
        rst = .false.
        print "(A)", "TEST FAILED: test_multiple_impute -1"
    end if

    ! Test 2 - the observed values must be left untouched
    do k = 1, nk
        do i = 1, n
            do j = 1, m
                if (is_missing(x(i,j))) cycle
                if (.not.assert(x(i,j), xc(i,j,k), 1.0d-10)) then
                    rst = .false.
                    print "(A)", "TEST FAILED: test_multiple_impute -2"
                    exit
                end if
            end do
        end do
    end do

    ! Test 3 - the imputations must differ from one another; the whole point of
    ! multiple imputation is to capture the uncertainty in the filled values
    varies = .false.
    do i = 1, n
        if (.not.is_missing(x(i,2))) cycle
        if (abs(xc(i,2,1) - xc(i,2,nk)) > 1.0d-10) then
            varies = .true.
            exit
        end if
    end do
    if (.not.varies) then
        rst = .false.
        print "(A)", "TEST FAILED: test_multiple_impute -3"
    end if

    ! Test 4 - each completed data set must recover the population means
    do k = 1, nk
        do j = 1, m
            avg(j) = sum(xc(:,j,k)) / real(n, real64)
        end do
        if (.not.assert(mtrue, avg, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_multiple_impute -4"
            exit
        end if
    end do
end function

! ------------------------------------------------------------------------------
function test_pool_imputations() result(rst)
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-10
    real(real64), parameter :: est(4) = [1.0d0, 2.0d0, 3.0d0, 4.0d0]
    real(real64), parameter :: vars(4) = [5.0d-1, 5.0d-1, 5.0d-1, 5.0d-1]

    ! Local Variables
    real(real64) :: qbar, tvar, dof, b, r, ans

    ! Initialization
    rst = .true.
    call pool_imputations(est, vars, qbar, tvar, dof)

    ! Test 1 - the pooled estimate is the mean of the estimates
    if (.not.assert(2.5d0, qbar, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_pool_imputations -1"
    end if

    ! Test 2 - the total variance is the within-imputation variance plus the
    ! scaled between-imputation variance
    b = sum((est - 2.5d0)**2) / 3.0d0
    ans = 5.0d-1 + 1.25d0 * b
    if (.not.assert(ans, tvar, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_pool_imputations -2"
    end if

    ! Test 3 - the degrees of freedom
    r = 1.25d0 * b / 5.0d-1
    ans = 3.0d0 * (1.0d0 + 1.0d0 / r)**2
    if (.not.assert(ans, dof, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_pool_imputations -3"
    end if

    ! Test 4 - a single imputation carries no between-imputation variance
    call pool_imputations(est(1:1), vars(1:1), qbar, tvar)
    if (.not.assert(est(1), qbar, tol) .or. .not.assert(vars(1), tvar, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_pool_imputations -4"
    end if
end function

! ------------------------------------------------------------------------------
end module

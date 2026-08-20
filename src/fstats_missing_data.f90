module fstats_missing_data
    !! Provides routines for dealing with data sets containing missing 
    !! observations.  Missing values are denoted by IEEE quiet NaN's.
    !!
    !! Three imputation strategies are provided.
    !!
    !! - k-nearest neighbors imputation ([[knn_impute]]): a nonparametric, 
    !!   donor-based approach.
    !!
    !! - Likelihood-based imputation ([[em_impute]]): maximum-likelihood 
    !!   estimation of the mean vector and covariance matrix of a multivariate
    !!   normal population by means of the expectation-maximization (EM) 
    !!   algorithm.  The missing values are replaced by their conditional
    !!   expectations.
    !!
    !! - Multiple imputation ([[multiple_impute]]): Bayesian data augmentation
    !!   under a multivariate normal model that generates several completed
    !!   copies of the data set such that the uncertainty introduced by the
    !!   imputation process can be propagated through subsequent analyses.  The
    !!   results of the analyses of each completed data set can be combined by
    !!   means of Rubin's rules ([[pool_imputations]]).
    use iso_fortran_env
    use ieee_arithmetic, only : ieee_value, ieee_is_nan, ieee_quiet_nan
    use fstats_errors
    use fstats_sampling, only : box_muller_sample
    use linalg, only : cholesky_factor, solve_cholesky
    implicit none
    private
    public :: missing_value
    public :: is_missing
    public :: knn_impute
    public :: em_impute
    public :: multiple_impute
    public :: pool_imputations

contains
! ******************************************************************************
! UTILITIES
! ------------------------------------------------------------------------------
pure function missing_value() result(rst)
    !! Returns the value used to flag a missing observation.  The value is an
    !! IEEE quiet NaN.
    real(real64) :: rst
        !! The missing-value flag.

    rst = ieee_value(0.0d0, ieee_quiet_nan)
end function

! ------------------------------------------------------------------------------
elemental function is_missing(x) result(rst)
    !! Determines if the supplied value represents a missing observation.
    real(real64), intent(in) :: x
        !! The value to test.
    logical :: rst
        !! Returns true if the value is missing; else, false.

    rst = ieee_is_nan(x)
end function

! ******************************************************************************
! K-NEAREST NEIGHBORS IMPUTATION
! ------------------------------------------------------------------------------
subroutine knn_impute(x, xc, k, weighted)
    !! Imputes missing values by means of a k-nearest neighbors approach.
    !!
    !! For each observation containing a missing value, the distance to every
    !! other observation is computed using only those variables observed in
    !! both rows.  The distance is a standardized, root-mean-square difference
    !! such that rows sharing differing numbers of observed variables remain
    !! comparable.  Each missing entry is then replaced by the mean of the
    !! corresponding entries of the k nearest rows for which that variable is
    !! observed.  If no suitable donor exists, the variable mean, computed from
    !! the observed values, is used.
    real(real64), intent(in), dimension(:,:) :: x
        !! An N-by-M matrix containing N observations of M variables.  Missing
        !! entries must be denoted by NaN's (see [[missing_value]]).
    real(real64), intent(out), dimension(:,:) :: xc
        !! An N-by-M matrix where the completed data set will be written.
    integer(int32), intent(in), optional :: k
        !! An optional input specifying the number of neighbors to use.  The
        !! default is 5.
    logical, intent(in), optional :: weighted
        !! An optional input that, if set to true, weights each donor by the
        !! inverse of its distance.  The default is true.

    ! Local Variables
    logical :: wgt
    logical, allocatable, dimension(:) :: used
    logical, allocatable, dimension(:,:) :: miss
    integer(int32) :: i, j, p, q, n, m, nk, nu, id, cnt
    real(real64) :: ds, dmin, num, den, w, big, eps
    real(real64), allocatable, dimension(:) :: mu, sd, dist

    ! Initialization
    n = size(x, 1)
    m = size(x, 2)
    big = huge(0.0d0)
    eps = sqrt(epsilon(0.0d0))
    if (present(k)) then
        nk = k
    else
        nk = 5
    end if
    if (present(weighted)) then
        wgt = weighted
    else
        wgt = .true.
    end if

    ! Input Checking
    if (size(xc, 1) /= n .or. size(xc, 2) /= m) error stop FS_MATRIX_SIZE_ERROR
    if (nk < 1) error stop FS_INVALID_INPUT_ERROR

    ! Compute the scaling factors from the observed data
    allocate(miss(n, m), used(n), dist(n), mu(m), sd(m))
    miss = is_missing(x)
    do j = 1, m
        if (all(miss(:,j))) error stop FS_INVALID_INPUT_ERROR
    end do
    call column_stats(x, miss, mu, sd)
    sd = sqrt(sd)
    where (sd <= 0.0d0) sd = 1.0d0

    ! Process
    do i = 1, n
        if (.not.any(miss(i,:))) then
            xc(i,:) = x(i,:)
            cycle
        end if

        ! Compute the distance from this row to every other row
        do p = 1, n
            dist(p) = big
            if (p == i) cycle
            ds = 0.0d0
            nu = 0
            do j = 1, m
                if (miss(i,j) .or. miss(p,j)) cycle
                ds = ds + ((x(i,j) - x(p,j)) / sd(j))**2
                nu = nu + 1
            end do
            if (nu > 0) dist(p) = sqrt(ds / real(nu, real64))
        end do

        ! Fill in each missing entry
        do j = 1, m
            if (.not.miss(i,j)) then
                xc(i,j) = x(i,j)
                cycle
            end if
            used = .false.
            num = 0.0d0
            den = 0.0d0
            cnt = 0
            do q = 1, nk
                ! Locate the nearest, as yet unused, row for which this
                ! variable is observed
                id = 0
                dmin = big
                do p = 1, n
                    if (used(p) .or. miss(p,j)) cycle
                    if (dist(p) < dmin) then
                        dmin = dist(p)
                        id = p
                    end if
                end do
                if (id == 0) exit
                used(id) = .true.
                if (wgt) then
                    w = 1.0d0 / (dmin + eps)
                else
                    w = 1.0d0
                end if
                num = num + w * x(id,j)
                den = den + w
                cnt = cnt + 1
            end do
            if (cnt == 0) then
                xc(i,j) = mu(j)
            else
                xc(i,j) = num / den
            end if
        end do
    end do
end subroutine

! ******************************************************************************
! LIKELIHOOD-BASED (EM) IMPUTATION
! ------------------------------------------------------------------------------
subroutine em_impute(x, xc, mu, sigma, maxiter, tol, niter)
    !! Employs the expectation-maximization (EM) algorithm to compute the
    !! maximum-likelihood estimates of the mean vector and covariance matrix of
    !! a multivariate normal population from a data set containing missing
    !! values.  The missing values are then replaced by their conditional
    !! expectations given the observed values and the converged parameter
    !! estimates.
    !!
    !! The estimator assumes the data are missing at random (MAR).
    real(real64), intent(in), dimension(:,:) :: x
        !! An N-by-M matrix containing N observations of M variables.  Missing
        !! entries must be denoted by NaN's (see [[missing_value]]).
    real(real64), intent(out), dimension(:,:) :: xc
        !! An N-by-M matrix where the completed data set will be written.
    real(real64), intent(out), optional, dimension(:) :: mu
        !! An optional M-element array where the maximum-likelihood estimate of
        !! the mean vector will be written.
    real(real64), intent(out), optional, dimension(:,:) :: sigma
        !! An optional M-by-M matrix where the maximum-likelihood estimate of
        !! the covariance matrix will be written.
    integer(int32), intent(in), optional :: maxiter
        !! An optional input specifying the maximum number of iterations to
        !! allow.  The default is 500.
    real(real64), intent(in), optional :: tol
        !! An optional input specifying the convergence tolerance.  The
        !! iteration terminates once the largest change in any parameter falls
        !! below this value.  The default is 1.0e-8.
    integer(int32), intent(out), optional :: niter
        !! An optional output containing the number of iterations performed.

    ! Local Variables
    logical, allocatable, dimension(:,:) :: miss
    integer(int32) :: i, n, m, iter, mxit
    real(real64) :: delta, tolerance
    real(real64), allocatable, dimension(:) :: mv, mnew, xf
    real(real64), allocatable, dimension(:,:) :: sv, snew, cv, t2

    ! Initialization
    n = size(x, 1)
    m = size(x, 2)
    if (present(maxiter)) then
        mxit = maxiter
    else
        mxit = 500
    end if
    if (present(tol)) then
        tolerance = tol
    else
        tolerance = 1.0d-8
    end if

    ! Input Checking
    if (size(xc, 1) /= n .or. size(xc, 2) /= m) error stop FS_MATRIX_SIZE_ERROR
    if (present(mu)) then
        if (size(mu) /= m) error stop FS_ARRAY_SIZE_ERROR
    end if
    if (present(sigma)) then
        if (size(sigma, 1) /= m .or. size(sigma, 2) /= m) then
            error stop FS_MATRIX_SIZE_ERROR
        end if
    end if
    if (mxit < 1) error stop FS_INVALID_INPUT_ERROR
    if (n < 2) error stop FS_UNDERDEFINED_PROBLEM_ERROR

    ! Establish the starting estimates
    allocate(miss(n, m), mv(m), mnew(m), xf(m), sv(m, m), snew(m, m), &
        cv(m, m), t2(m, m))
    miss = is_missing(x)
    call starting_estimates(x, miss, mv, sv)

    ! The EM iteration
    do iter = 1, mxit
        ! E-step: accumulate the expected sufficient statistics
        mnew = 0.0d0
        t2 = 0.0d0
        do i = 1, n
            call conditional_step(x(i,:), miss(i,:), mv, sv, .false., xf, cv)
            mnew = mnew + xf
            t2 = t2 + outer_product(xf, xf) + cv
        end do

        ! M-step
        mnew = mnew / real(n, real64)
        snew = t2 / real(n, real64) - outer_product(mnew, mnew)
        snew = 0.5d0 * (snew + transpose(snew))

        delta = max(maxval(abs(mnew - mv)), maxval(abs(snew - sv)))
        mv = mnew
        sv = snew
        if (delta < tolerance) exit
    end do
    iter = min(iter, mxit)

    ! Replace the missing values with their conditional expectations
    do i = 1, n
        call conditional_step(x(i,:), miss(i,:), mv, sv, .false., xf, cv)
        xc(i,:) = xf
    end do

    ! Output
    if (present(mu)) mu = mv
    if (present(sigma)) sigma = sv
    if (present(niter)) niter = iter
end subroutine

! ******************************************************************************
! MULTIPLE IMPUTATION
! ------------------------------------------------------------------------------
subroutine multiple_impute(x, xc, nburn, maxiter, tol)
    !! Generates multiple completed copies of a data set containing missing
    !! values by means of Bayesian data augmentation under a multivariate
    !! normal model with the noninformative Jeffreys prior.
    !!
    !! The EM algorithm ([[em_impute]]) supplies the starting parameter
    !! estimates.  Thereafter, each imputation is generated by alternating an
    !! imputation step, wherein the missing values are drawn from their
    !! conditional distribution given the observed values and the current
    !! parameters, with a posterior step, wherein new parameters are drawn from
    !! their posterior distribution given the completed data set.
    !!
    !! The analyses of the resulting data sets may be combined by means of
    !! [[pool_imputations]].
    real(real64), intent(in), dimension(:,:) :: x
        !! An N-by-M matrix containing N observations of M variables.  Missing
        !! entries must be denoted by NaN's (see [[missing_value]]).
    real(real64), intent(out), dimension(:,:,:) :: xc
        !! An N-by-M-by-K array where the K completed data sets will be
        !! written.
    integer(int32), intent(in), optional :: nburn
        !! An optional input specifying the number of data augmentation cycles
        !! to perform between successive imputations.  The default is 20.
    integer(int32), intent(in), optional :: maxiter
        !! An optional input specifying the maximum number of EM iterations to
        !! allow when computing the starting estimates.  The default is 500.
    real(real64), intent(in), optional :: tol
        !! An optional input specifying the convergence tolerance used when
        !! computing the EM starting estimates.  The default is 1.0e-8.

    ! Local Variables
    logical, allocatable, dimension(:,:) :: miss
    integer(int32) :: i, j, n, m, nk, imp, it, nb, df
    real(real64), allocatable, dimension(:) :: mv, xbar, xf
    real(real64), allocatable, dimension(:,:) :: sv, cv, xfull, s, lc

    ! Initialization
    n = size(x, 1)
    m = size(x, 2)
    nk = size(xc, 3)
    if (present(nburn)) then
        nb = nburn
    else
        nb = 20
    end if

    ! Input Checking
    if (size(xc, 1) /= n .or. size(xc, 2) /= m) error stop FS_MATRIX_SIZE_ERROR
    if (nk < 1 .or. nb < 1) error stop FS_INVALID_INPUT_ERROR

    ! The Wishart posterior requires more observations than variables
    df = n - 1
    if (df < m) error stop FS_UNDERDEFINED_PROBLEM_ERROR

    ! Compute the starting values by way of the EM algorithm
    allocate(miss(n, m), mv(m), xbar(m), xf(m), sv(m, m), cv(m, m), &
        xfull(n, m), s(m, m), lc(m, m))
    miss = is_missing(x)
    call em_impute(x, xfull, mu = mv, sigma = sv, maxiter = maxiter, tol = tol)

    ! Data augmentation
    do imp = 1, nk
        do it = 1, nb
            ! I-step: draw the missing values from their conditional
            ! distribution given the observed values and the current parameters
            do i = 1, n
                call conditional_step(x(i,:), miss(i,:), mv, sv, .true., xf, cv)
                xfull(i,:) = xf
            end do

            ! P-step: draw new parameters from their posterior distribution
            do j = 1, m
                xbar(j) = sum(xfull(:,j)) / real(n, real64)
            end do
            s = 0.0d0
            do i = 1, n
                s = s + outer_product(xfull(i,:) - xbar, xfull(i,:) - xbar)
            end do
            call draw_inverse_wishart(s, df, sv)
            call spd_cholesky(sv, lc)
            mv = xbar + matmul(lc, std_normal(m)) / sqrt(real(n, real64))
        end do
        xc(:,:,imp) = xfull
    end do
end subroutine

! ------------------------------------------------------------------------------
pure subroutine pool_imputations(est, vars, qbar, tvar, dof)
    !! Combines the results of the analyses of multiply-imputed data sets by
    !! means of Rubin's rules.
    real(real64), intent(in), dimension(:) :: est
        !! A K-element array containing the parameter estimate obtained from
        !! each of the K imputed data sets.
    real(real64), intent(in), dimension(:) :: vars
        !! A K-element array containing the estimated variance of the parameter
        !! obtained from each of the K imputed data sets.
    real(real64), intent(out) :: qbar
        !! The pooled parameter estimate.
    real(real64), intent(out) :: tvar
        !! The total variance of the pooled estimate.  The total variance is
        !! the sum of the within-imputation variance and the appropriately
        !! scaled between-imputation variance.
    real(real64), intent(out), optional :: dof
        !! An optional output containing the degrees of freedom associated with
        !! the pooled estimate.

    ! Local Variables
    integer(int32) :: k
    real(real64) :: ubar, b, r, kr

    ! Initialization
    k = size(est)
    kr = real(k, real64)

    ! Quick Return
    if (k < 2) then
        if (k == 1) then
            qbar = est(1)
            tvar = vars(1)
        else
            qbar = 0.0d0
            tvar = 0.0d0
        end if
        if (present(dof)) dof = 0.0d0
        return
    end if

    ! Process
    qbar = sum(est) / kr
    ubar = sum(vars) / kr
    b = sum((est - qbar)**2) / (kr - 1.0d0)
    tvar = ubar + (1.0d0 + 1.0d0 / kr) * b
    if (present(dof)) then
        if (b <= 0.0d0 .or. ubar <= 0.0d0) then
            dof = huge(0.0d0)
        else
            r = (1.0d0 + 1.0d0 / kr) * b / ubar
            dof = (kr - 1.0d0) * (1.0d0 + 1.0d0 / r)**2
        end if
    end if
end subroutine

! ******************************************************************************
! PRIVATE SUPPORTING ROUTINES
! ------------------------------------------------------------------------------
pure subroutine column_stats(x, miss, mu, var)
    !! Computes the mean and variance of each column using only the observed
    !! values.
    real(real64), intent(in), dimension(:,:) :: x
    logical, intent(in), dimension(:,:) :: miss
    real(real64), intent(out), dimension(:) :: mu
    real(real64), intent(out), dimension(:) :: var

    ! Local Variables
    integer(int32) :: i, j, n, m, cnt

    n = size(x, 1)
    m = size(x, 2)
    do j = 1, m
        mu(j) = 0.0d0
        var(j) = 0.0d0
        cnt = 0
        do i = 1, n
            if (miss(i,j)) cycle
            mu(j) = mu(j) + x(i,j)
            cnt = cnt + 1
        end do
        if (cnt == 0) cycle
        mu(j) = mu(j) / real(cnt, real64)
        if (cnt < 2) cycle
        do i = 1, n
            if (miss(i,j)) cycle
            var(j) = var(j) + (x(i,j) - mu(j))**2
        end do
        var(j) = var(j) / real(cnt - 1, real64)
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine starting_estimates(x, miss, mu, sigma)
    !! Computes the starting values for the EM iteration.  The mean vector is
    !! the vector of available-case column means, and the covariance matrix is
    !! taken as diagonal to guarantee a positive-definite starting point.
    real(real64), intent(in), dimension(:,:) :: x
    logical, intent(in), dimension(:,:) :: miss
    real(real64), intent(out), dimension(:) :: mu
    real(real64), intent(out), dimension(:,:) :: sigma

    ! Local Variables
    integer(int32) :: j, m
    real(real64), allocatable, dimension(:) :: var

    ! A variable with no observed values cannot be imputed
    m = size(x, 2)
    do j = 1, m
        if (all(miss(:,j))) error stop FS_INVALID_INPUT_ERROR
    end do

    allocate(var(m))
    call column_stats(x, miss, mu, var)
    sigma = 0.0d0
    do j = 1, m
        if (var(j) <= 0.0d0) then
            sigma(j,j) = 1.0d0
        else
            sigma(j,j) = var(j)
        end if
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine conditional_step(xrow, miss, mu, sigma, draw, xf, cvar)
    !! Computes the conditional mean and conditional covariance of the missing
    !! entries of a single observation given its observed entries.  If
    !! requested, the missing entries are sampled at random from the resulting
    !! conditional normal distribution rather than being set to the conditional
    !! mean.
    real(real64), intent(in), dimension(:) :: xrow
    logical, intent(in), dimension(:) :: miss
    real(real64), intent(in), dimension(:) :: mu
    real(real64), intent(in), dimension(:,:) :: sigma
    logical, intent(in) :: draw
    real(real64), intent(out), dimension(:) :: xf
        !! The completed observation.
    real(real64), intent(out), dimension(:,:) :: cvar
        !! The contribution to the expected sufficient statistics arising from
        !! the conditional variance.  All entries are zero except those on the
        !! rows and columns corresponding to the missing entries.  The matrix
        !! is zeroed if DRAW is true.

    ! Local Variables
    integer(int32) :: i, m, nm, no, io, im
    integer(int32), allocatable, dimension(:) :: iobs, imis
    real(real64), allocatable, dimension(:) :: xhat
    real(real64), allocatable, dimension(:,:) :: soo, som, smm, z, lc, cmm

    ! Initialization
    m = size(xrow)
    nm = count(miss)
    no = m - nm
    xf = xrow
    cvar = 0.0d0

    ! Quick Return - nothing is missing
    if (nm == 0) return

    ! Handle the case where the entire observation is missing
    if (no == 0) then
        if (draw) then
            allocate(lc(m, m))
            call spd_cholesky(sigma, lc)
            xf = mu + matmul(lc, std_normal(m))
        else
            xf = mu
            cvar = sigma
        end if
        return
    end if

    ! Partition the observation
    allocate(iobs(no), imis(nm), xhat(nm), soo(no, no), som(no, nm), &
        smm(nm, nm), z(no, nm), cmm(nm, nm), lc(no, no))
    io = 0
    im = 0
    do i = 1, m
        if (miss(i)) then
            im = im + 1
            imis(im) = i
        else
            io = io + 1
            iobs(io) = i
        end if
    end do
    soo = sigma(iobs,iobs)
    som = sigma(iobs,imis)
    smm = sigma(imis,imis)

    ! The conditional mean is mu(m) + S(m,o) * inv(S(o,o)) * (x(o) - mu(o)),
    ! and the conditional covariance is S(m,m) - S(m,o) * inv(S(o,o)) * S(o,m).
    ! Both are formed from Z = inv(S(o,o)) * S(o,m).
    call spd_cholesky(soo, lc)
    z = solve_cholesky(.false., lc, som)
    xhat = mu(imis) + matmul(transpose(z), xrow(iobs) - mu(iobs))
    cmm = smm - matmul(transpose(som), z)
    cmm = 0.5d0 * (cmm + transpose(cmm))

    ! Store the results
    if (draw) then
        deallocate(lc)
        allocate(lc(nm, nm))
        call spd_cholesky(cmm, lc)
        xf(imis) = xhat + matmul(lc, std_normal(nm))
    else
        xf(imis) = xhat
        cvar(imis,imis) = cmm
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine draw_inverse_wishart(s, df, rst)
    !! Draws a matrix at random from an inverse Wishart distribution with scale
    !! matrix S and DF degrees of freedom by means of the Bartlett 
    !! decomposition.
    real(real64), intent(in), dimension(:,:) :: s
    integer(int32), intent(in) :: df
    real(real64), intent(out), dimension(:,:) :: rst

    ! Local Variables
    integer(int32) :: i, m
    real(real64), allocatable, dimension(:) :: z
    real(real64), allocatable, dimension(:,:) :: sinv, c, a, w

    ! A draw W from a Wishart distribution with DF degrees of freedom and scale
    ! matrix inv(S) is formed as W = (C * A) * (C * A)**T, where C is the
    ! Cholesky factor of inv(S) and A is the Bartlett factor.  The desired
    ! inverse Wishart draw is then inv(W).
    m = size(s, 1)
    allocate(a(m, m), c(m, m), sinv(m, m), w(m, m))
    sinv = spd_inverse(s)
    call spd_cholesky(sinv, c)

    a = 0.0d0
    do i = 1, m
        z = std_normal(df - i + 1)
        a(i,i) = sqrt(sum(z**2))
    end do
    do i = 2, m
        a(i,1:i-1) = std_normal(i - 1)
    end do

    ! W = (C * A) * (C * A)**T is already in factored form
    w = matmul(c, a)
    rst = solve_cholesky(.false., w, identity_matrix(m))
    rst = 0.5d0 * (rst + transpose(rst))
end subroutine

! ------------------------------------------------------------------------------
function std_normal(n) result(rst)
    !! Generates N standard normal deviates.
    integer(int32), intent(in) :: n
    real(real64), allocatable, dimension(:) :: rst

    rst = box_muller_sample(0.0d0, 1.0d0, n)
end function

! ------------------------------------------------------------------------------
pure function outer_product(a, b) result(rst)
    !! Computes the outer product of two vectors.
    real(real64), intent(in), dimension(:) :: a, b
    real(real64), allocatable, dimension(:,:) :: rst

    ! Local Variables
    integer(int32) :: i, j

    allocate(rst(size(a), size(b)))
    do j = 1, size(b)
        do i = 1, size(a)
            rst(i,j) = a(i) * b(j)
        end do
    end do
end function

! ------------------------------------------------------------------------------
pure function identity_matrix(n) result(rst)
    !! Constructs an N-by-N identity matrix.
    integer(int32), intent(in) :: n
    real(real64), allocatable, dimension(:,:) :: rst

    ! Local Variables
    integer(int32) :: i

    allocate(rst(n, n), source = 0.0d0)
    do i = 1, n
        rst(i,i) = 1.0d0
    end do
end function

! ------------------------------------------------------------------------------
subroutine spd_cholesky(a, l)
    !! Computes the lower-triangular Cholesky factorization of a symmetric
    !! matrix.  A small multiple of the identity matrix is added to the
    !! diagonal if the matrix is not sufficiently positive-definite.
    real(real64), intent(in), dimension(:,:) :: a
    real(real64), intent(out), dimension(:,:) :: l

    ! Local Variables
    integer(int32) :: i, n, attempt
    real(real64) :: jitter, scaling
    real(real64), allocatable, dimension(:,:) :: b

    n = size(a, 1)
    b = a
    scaling = 0.0d0
    do i = 1, n
        scaling = scaling + abs(a(i,i))
    end do
    scaling = max(scaling / real(n, real64), 1.0d-16)
    jitter = 0.0d0
    do attempt = 1, 10
        l = cholesky_factor(b, .false.)
        if (is_valid_factor(l)) return
        if (jitter == 0.0d0) then
            jitter = 1.0d-12 * scaling
        else
            jitter = 1.0d2 * jitter
        end if
        b = a
        do i = 1, n
            b(i,i) = b(i,i) + jitter
        end do
    end do

    ! The matrix could not be factored
    error stop FS_INVALID_INPUT_ERROR
end subroutine

! ------------------------------------------------------------------------------
pure function is_valid_factor(l) result(rst)
    !! Determines if a Cholesky factorization succeeded.  A nonpositive or NaN
    !! diagonal term denotes a matrix that is not positive-definite.
    real(real64), intent(in), dimension(:,:) :: l
    logical :: rst

    ! Local Variables
    integer(int32) :: i

    rst = .true.
    do i = 1, size(l, 1)
        if (.not.(l(i,i) > 0.0d0)) then
            rst = .false.
            return
        end if
    end do
end function

! ------------------------------------------------------------------------------
function spd_inverse(a) result(rst)
    !! Computes the inverse of a symmetric, positive-definite matrix.
    real(real64), intent(in), dimension(:,:) :: a
    real(real64), allocatable, dimension(:,:) :: rst

    ! Local Variables
    integer(int32) :: n
    real(real64), allocatable, dimension(:,:) :: l

    n = size(a, 1)
    allocate(l(n, n))
    call spd_cholesky(a, l)
    rst = solve_cholesky(.false., l, identity_matrix(n))
end function

! ------------------------------------------------------------------------------
end module
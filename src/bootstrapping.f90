submodule (fstats) bootstrapping
    implicit none

! REFERENCES:
! - https://medium.com/@m21413108/bootstrapping-maximum-entropy-non-parametric-boot-python-3b1e23ea589d
! - https://cran.r-project.org/web/packages/meboot/vignettes/meboot.pdf
! - https://gist.github.com/christianjauregui/314456688a3c2fead43a48be3a47dad6

contains
! ------------------------------------------------------------------------------
module subroutine bs_linear_least_squares_real64(order, intercept, x, y, &
    coeffs, ymod, resid, nsamples, stats, bias, alpha, err)
    ! Arguments
    integer(int32), intent(in) :: order
    logical, intent(in) :: intercept
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: coeffs
    real(real64), intent(out), dimension(:) :: ymod
    real(real64), intent(out), dimension(:) :: resid
    integer(int32), intent(in), optional :: nsamples
    type(bootstrap_regression_statistics), intent(out), optional, &
        dimension(:) :: stats
    real(real64), intent(out), optional, dimension(:) :: bias
    real(real64), intent(in), optional :: alpha
    class(errors), intent(inout), optional, target :: err

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: half = 0.5d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, j, i1, i2, n, ns, nc, ncoeffs, flag
    real(real64) :: eps, alph, ms
    real(real64), allocatable, dimension(:,:) :: allcoeffs, allresid, allf, &
        ally
    type(t_distribution) :: dist
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
        alph = 0.05d0
    end if
    n = size(x)
    ncoeffs = order + 1
    nc = order
    i1 = floor(alph * ns, int32)
    i2 = ns - i1 + 1
    if (intercept) nc = nc + 1
    dist%dof = real(ns - nc)

    ! TO DO:
    ! - Input Checking
    ! - Verify the size of BIAS & ensure it's the appropriate size along with
    !   COEFF

    ! Compute the fit
    call linear_least_squares(order, intercept, x, y, coeffs, &
        ymod, resid, alpha = alph, err = errmgr)
    if (errmgr%has_error_occurred()) return
    
    ! Memory Allocations
    allocate(allcoeffs(ncoeffs, ns), allresid(n, ns), allf(n, ns), ally(n, ns - 1), &
        stat = flag)
    if (flag /= 0) then
        ! TO DO: Memory Error
    end if
    allcoeffs(:,1) = coeffs

    ! Establish the epsilon term - to do, come up with a better resampling approach
    eps = standard_deviation(y) / sqrt(real(n))
    call random_init(.false., .true.)
    call random_number(ally)
    ally = eps * (ally - half)

    ! Cycle over each data set
    do i = 2, ns
        ! Compute the fit of the perturbed data set
        ally(:,i-1) = ally(:,i-1) + y
        call linear_least_squares(order, intercept, x, ally(:,i-1), &
            allcoeffs(:,i), allf(:,i), allresid(:,i), alpha = alph)
    end do
    
    ! Perform statistics calculations, if needed
    if (present(stats)) then
        ! Update the relevant statistical metrics for each coefficient based
        ! upon the actual distribution
        j = 1
        if (intercept) j = 0
        do i = 1, nc
            j = j + 1
            ms = trimmed_mean(allcoeffs(j,:), p = half * alph)
            ! As we have a distribution of mean values, the standard deviation
            ! of this population yields the standard error estimate for the
            ! overall problem
            stats(i)%standard_error = standard_deviation(allcoeffs(j,:))
            ! As before, this is a distribution of mean values.  The CI can
            ! be directly estimated by considering the values of the bottom
            ! alpha/2 and top alpha/2 terms.
            stats(i)%upper_confidence_interval = allcoeffs(j,i2)
            stats(i)%lower_confidence_interval = allcoeffs(j,i1)
            stats(i)%t_statistic = coeffs(i) / stats(i)%standard_error
            stats(i)%probability = regularized_beta(half * dist%dof, half, &
                dist%dof / (dist%dof + (stats(i)%t_statistic)**2))
        end do
    end if

    ! Compute the bias for each parameter, if needed
    if (present(bias)) then
        ! Verify the size of the array

        ! Perform the calculations
        do i = 1, nc
            bias(i) = coeffs(i) - mean(allcoeffs(i,:))
        end do
    end if
end subroutine

! ------------------------------------------------------------------------------
end submodule
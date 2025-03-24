module fstats_sampling
    use iso_fortran_env
    use linalg, only : sort
    use fstats_distributions
    implicit none
    private
    public :: box_muller_sample
    public :: rejection_sample
    public :: sample_normal_multivariate

    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: twopi = 2.0d0 * pi

    ! interface box_muller_sample
    !     !! Generates random, normally distributed values via the Box-Muller 
    !     !! transform.
    !     module procedure :: box_muller_sample_scalar
    !     module procedure :: box_muller_array
    ! end interface
contains
! ------------------------------------------------------------------------------
function box_muller_sample_driver(mu, sigma, prev) result(rst)
    !! Driver routine for computing a single Box-Muller sample pair.
    real(real64), intent(in) :: mu
    !! The mean of the distribution.    
    real(real64), intent(in) :: sigma
        !! The standard deviation of the distribution.
    real(real64), intent(inout) :: prev
        !! The extra value from the previous set of samples.  On output, this
        !! value is updated.
    real(real64) :: rst
        !! The resulting sample.

    ! Local Variables
    real(real64) :: u1, u2, v1, v2, rsq, fac

    ! Process
    if (prev == 0.0d0) then
        rsq = 1.0d0
        do while (rsq >= 1.0d0 .or.rsq == 0.0d0)
            call random_number(u1)
            call random_number(u2)
            v1 = 2.0d0 * u1 - 1.0d0
            v2 = 2.0d0 * u2 - 1.020
            rsq = v1**2 + v2**2
        end do
        fac = sqrt(-2.0d0 * log(rsq) / rsq)
        prev = v1 * fac
        rst = mu + sigma * v2 * fac
    else
        fac = prev
        rst = mu + sigma * fac
        prev = 0.0d0
    end if
end function

! ------------------------------------------------------------------------------
function box_muller_sample(mu, sigma, n) result(rst)
    !! Utilizes the Box-Muller transformation approach to generate the requested
    !! number of random samples from a normal distribution of the specified mean
    !! and standard deviation.
    real(real64), intent(in) :: mu
        !! The mean of the distribution.
    real(real64), intent(in) :: sigma
        !! The standard deviation of the distribution.
    integer(int32), intent(in) :: n
        !! The number of random samples to generate.
    real(real64), allocatable, dimension(:) :: rst
        !! An N-element array containing the values.

    ! Local Variables
    integer(int32) :: i
    real(real64) :: prev

    ! Process
    if (n < 1) then
        allocate(rst(0))
        return
    end if
    allocate(rst(n))
    prev = 0.0d0
    do i = 1, n
        rst(i) = box_muller_sample_driver(mu, sigma, prev)
    end do
end function

! ******************************************************************************
! REJECTION SAMPLING
! ------------------------------------------------------------------------------
function rejection_sample(tdist, n, xmin, xmax) result(rst)
    !! Uses rejection sampling to randomly sample a target distribution.
    class(distribution), intent(in) :: tdist
        !! The distribution to sample
    integer(int32), intent(in) :: n
        !! The number of samples to make.
    real(real64), intent(in) :: xmin
        !! The minimum range to explore.
    real(real64), intent(in) :: xmax
        !! The maximum range to explore.
    real(real64), allocatable, dimension(:) :: rst
        !! An N-element array containing the N samples from the 
        !! distribution.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: c_start = 1.01d0

    ! Local Variables
    integer(int32) :: i, j, jmax
    real(real64) :: u, c, g, f, rng

    ! Quick Return
    if (n < 1) then
        allocate(rst(0), source = zero)
    end if

    ! Process
    i = 0
    j = 0
    jmax = min(1000 * n, huge(j))  ! Guard against insanity
    rng = xmax - xmin
    c = c_start
    allocate(rst(n), source = zero)
    do while (i < n)
        ! Update the acceptance threshold
        call random_number(u)

        ! Sample from the proposal distribution
        call random_number(g)
        g = g * rng + xmin

        ! Sample the target distribution
        f = tdist%pdf(g)

        ! Test
        if (u <= f / (c * g)) then
            i = i + 1
            rst(i) = g
        end if

        ! Update C
        c = max(c, f / g)

        ! Update the infinite loop guard variable
        j = j + 1
        if (j == jmax) exit
    end do
end function

! ******************************************************************************
! MULTIVARIATE SAMPLING
! ------------------------------------------------------------------------------
function sample_normal_multivariate(dist) result(rst)
    !! Samples a multivariate normal distribution such that \(\vec{x} = 
    !! \vec{mu} + L \vec{u}\), where \(L\) is the lower form of the Cholesky 
    !! factorization of the covariance matrix, and \(\vec{u}\) is a randomly 
    !! generated vector that exists on the set \([-1 1]\)
    class(multivariate_normal_distribution), intent(in) :: dist
        !! The multivariate normal distribution to sample.
    real(real64), allocatable, dimension(:) :: rst
        !! The resulting vector.

    ! Local Variables
    integer(int32) :: n
    real(real64), allocatable, dimension(:) :: u
    real(real64), allocatable, dimension(:,:) :: L

    ! Initialization
    L = dist%get_cholesky_factored_matrix()
    n = size(L, 1)
    allocate(u(n))
    call random_number(u)   ! populating u from [0, 1].
    u = 2.0d0 * (u - 0.5d0) ! centering u around 0 over the range [-1, 1]

    ! Process
    rst = dist%get_means() + matmul(L, u)
end function

! ------------------------------------------------------------------------------
end module
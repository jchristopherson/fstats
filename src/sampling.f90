submodule (fstats) sampling
    use linalg, only : sort
    implicit none

    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: twopi = 2.0d0 * pi
    real(real64), parameter :: pi_f = 2.0 * acos(0.0)
    real(real64), parameter :: twopi_f = 2.0 * pi_f
contains
! ------------------------------------------------------------------------------
module function box_muller_sample_real64(mu, sigma) result(rst)
    ! Arguments
    real(real64), intent(in) :: mu
    real(real64), intent(in) :: sigma
    real(real64) :: rst(2)

    ! Parameters
    complex(real64), parameter :: j = (0.0d0, 1.0d0)

    ! Local Variables
    real(real64) :: u1, u2
    complex(real64) :: z

    ! Process
    call random_number(u1)
    call random_number(u2)
    z = sqrt(-log(u1)) * exp(j * twopi * u2)
    rst = [real(z, real64), aimag(z)]
end function

! ------------------------------------------------------------------------------
module function box_muller_array_real64(mu, sigma, n) result(rst)
    ! Arguments
    real(real64), intent(in) :: mu
    real(real64), intent(in) :: sigma
    integer(int32), intent(in) :: n
    real(real64), allocatable, dimension(:) :: rst

    ! Local Variables
    integer(int32) :: i

    ! Process
    if (n < 1) then
        allocate(rst(0))
        return
    end if
    allocate(rst(2 * n))
    do i = 1, n
        rst(2*i-1:2*i) = box_muller_sample(mu, sigma)
    end do
end function

! ******************************************************************************
! REJECTION SAMPLING
! ------------------------------------------------------------------------------
module function rejection_sample(tdist, n, xmin, xmax) result(rst)
    ! Arguments
    class(distribution), intent(in) :: tdist
    integer(int32), intent(in) :: n
    real(real64), intent(in) :: xmin, xmax
    real(real64), allocatable, dimension(:) :: rst

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
    do while (i <= n)
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

end submodule
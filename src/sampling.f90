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

! --------------------
module function box_muller_sample_real32(mu, sigma) result(rst)
    ! Arguments
    real(real32), intent(in) :: mu
    real(real32), intent(in) :: sigma
    real(real32) :: rst(2)

    ! Parameters
    complex(real32), parameter :: j = (0.0, 1.0)

    ! Local Variables
    real(real32) :: u1, u2
    complex(real32) :: z

    ! Process
    call random_number(u1)
    call random_number(u2)
    z = sqrt(-log(u1)) * exp(j * twopi_f * u2)
    rst = [real(z, real32), aimag(z)]
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

! --------------------
module function box_muller_array_real32(mu, sigma, n) result(rst)
    ! Arguments
    real(real32), intent(in) :: mu
    real(real32), intent(in) :: sigma
    integer(int32), intent(in) :: n
    real(real32), allocatable, dimension(:) :: rst

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
! MAXIMUM ENTROPY
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end submodule
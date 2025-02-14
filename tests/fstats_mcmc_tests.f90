module fstats_mcmc_tests
    use iso_fortran_env
    use fstats
    use fortran_test_helper
    implicit none
contains
! ------------------------------------------------------------------------------
function test_mh_push() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: nvars = 10
    integer(int32), parameter :: npts = 100

    ! Local Variables
    integer(int32) :: i
    real(real64) :: x(npts, nvars)
    real(real64), allocatable, dimension(:,:) :: chain
    type(metropolis_hastings) :: mcmc

    ! Initialization
    rst = .true.
    call random_number(x)

    ! Push items onto the stack
    do i = 1, npts
        call mcmc%push_new_state(x(i,:))
    end do

    ! Check the stack size
    if (mcmc%get_chain_length() /= npts) then
        rst = .false.
        print "(A)", "TEST FAILED: test_mh_push -1"
    end if

    ! Check the number of variables
    if (mcmc%get_state_variable_count() /= nvars) then
        rst = .false.
        print "(A)", "TEST FAILED: test_mh_push -2"
    end if

    ! Get the chain
    chain = mcmc%get_chain()
    if (.not.assert(x, chain)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_mh_push -3"
    end if
end function

! ------------------------------------------------------------------------------
subroutine exfun(x, p, f, stop)
    ! Arguments
    real(real64), intent(in) :: x(:), p(:)
    real(real64), intent(out) :: f(:)
    logical, intent(out) :: stop

    ! Function
    f = p(4) * x**3 + p(3) * x**2 + p(2) * x + p(1)

    ! No need to stop
    stop = .false.
end subroutine

! ------------------------------------------------------------------------------
function test_mh_proposal() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: mdl(4), xc(4), xp(4), x(21), y(21)
    type(metropolis_hastings) :: mcmc
    procedure(regression_function), pointer :: fcn

    ! Initialization
    rst = .true.
    fcn => exfun
    call random_number(mdl)
    call random_number(xc)
    
    x = [0.0d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.8d0, &
        0.9d0, 1.0d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0, 1.5d0, 1.6d0, 1.7d0, &
        1.8d0, 1.9d0, 2.0d0]
    y = [1.216737514d0, 1.250032542d0, 1.305579195d0, 1.040182335d0, &
        1.751867738d0, 1.109716707d0, 2.018141531d0, 1.992418729d0, &
        1.807916923d0, 2.078806005d0, 2.698801324d0, 2.644662712d0, &
        3.412756702d0, 4.406137221d0, 4.567156645d0, 4.999550779d0, &
        5.652854194d0, 6.784320119d0, 8.307936836d0, 8.395126494d0, &
        10.30252404d0]

    ! Initialize the MH object
    call mcmc%initialize(fcn, x, y, mdl)

    ! Generate a new proposal
    call mcmc%proposal(xc, xp)

    ! Ensure xc and xp are not equal
    if (assert(xc, xp)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_mh_proposal -1"
    end if
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
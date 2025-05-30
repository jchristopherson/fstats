module fstats_mcmc_tests
    use iso_fortran_env
    use fstats
    use fortran_test_helper
    implicit none

    type, extends(mcmc_target) :: test_mcmc_target
    contains
        procedure, public :: model => tmt_eval
    end type

contains
! ******************************************************************************
! TEST_MCMC_TARGET
! ------------------------------------------------------------------------------
subroutine tmt_eval(this, xdata, xc, y)
    class(test_mcmc_target), intent(in) :: this
    real(real64), intent(in), dimension(:) :: xdata
    real(real64), intent(in), dimension(:) :: xc
    real(real64), intent(out), dimension(:) :: y

    ! Linear Model: y = m * x + b
    y = xc(1) * xdata + xc(2)
end subroutine

! ******************************************************************************
! MCMC_TARGET TESTS
! ------------------------------------------------------------------------------
function test_mcmc_target_distributions() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: mean_ans_1, mean_ans_2, std_ans_1, std_ans_2
    type(normal_distribution) :: dist1, dist2
    type(test_mcmc_target) :: target
    class(distribution), pointer :: dist

    ! Initialization
    rst = .true.
    call random_number(mean_ans_1)
    call random_number(mean_ans_2)
    call random_number(std_ans_1)
    call random_number(std_ans_2)
    dist1%mean_value = mean_ans_1
    dist1%standard_deviation = std_ans_1
    dist2%mean_value = mean_ans_2
    dist2%standard_deviation = std_ans_2

    ! Define distributions for each model parameter
    call target%add_parameter(dist1)
    call target%add_parameter(dist2)

    ! Test 1
    if (target%get_parameter_count() /= 2) then
        rst = .false.
        print "(A)", "Test Failed: test_mcmc_target_distributions -1"
        return
    end if

    ! Test 2
    dist => target%get_parameter(1)
    if (.not.assert(dist%mean(), mean_ans_1)) then
        rst = .false.
        print "(A)", "Test Failed: test_mcmc_target_distributions -2"
    end if
    if (.not.assert(sqrt(dist%variance()), std_ans_1)) then
        rst = .false.
        print "(A)", "Test Failed: test_mcmc_target_distributions -3"
    end if

    ! Test 3
    dist => target%get_parameter(2)
    if (.not.assert(dist%mean(), mean_ans_2)) then
        rst = .false.
        print "(A)", "Test Failed: test_mcmc_target_distributions -4"
    end if
    if (.not.assert(sqrt(dist%variance()), std_ans_2)) then
        rst = .false.
        print "(A)", "Test Failed: test_mcmc_target_distributions -5"
    end if
end function

! ------------------------------------------------------------------------------
function test_mcmc_target_likelihood() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: ndata = 21
    real(real64), parameter :: var = 1.0d-2
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    integer(int32) :: i
    real(real64) :: v, prd, l, xdata(ndata), ydata(ndata), ymod(ndata), params(2)
    type(normal_distribution) :: dist, param1, param2
    type(test_mcmc_target) :: target

    ! Initialization
    rst = .true.
    xdata = [0.0d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.8d0, &
            0.9d0, 1.0d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0, 1.5d0, 1.6d0, 1.7d0, &
            1.8d0, 1.9d0, 2.0d0]
    ydata = [1.216737514d0, 1.250032542d0, 1.305579195d0, 1.040182335d0, &
            1.751867738d0, 1.109716707d0, 2.018141531d0, 1.992418729d0, &
            1.807916923d0, 2.078806005d0, 2.698801324d0, 2.644662712d0, &
            3.412756702d0, 4.406137221d0, 4.567156645d0, 4.999550779d0, &
            5.652854194d0, 6.784320119d0, 8.307936836d0, 8.395126494d0, &
            10.30252404d0]

    ! Use randomly generated model parameters
    call random_number(params)

    ! Set up the target
    param1%mean_value = params(1)
    param1%standard_deviation = sqrt(var)

    param2%mean_value = params(2)
    param2%standard_deviation = sqrt(var)

    call target%add_parameter(param1)
    call target%add_parameter(param2)
    
    ! Evaluate the model - store in ymod
    call target%model(xdata, params, ymod)

    ! Evaluate the distribution - this will be our answer term
    ! Use a log approach to avoid underflow/overflow issues
    dist%standard_deviation = sqrt(var)
    prd = 0.0d0
    do i = 1, ndata
        dist%mean_value = ymod(i)
        v = log(dist%pdf(ydata(i)))  ! evaluate the distribution at the measued y
        prd = prd + v
    end do
    prd = exp(prd)

    ! Compute the likelihood via the target type
    l = target%likelihood(xdata, ydata, params, var)

    ! Test
    if (.not.assert(prd, l, tol)) then
        rst = .false.
        print "(A)", "Test Failed: test_mcmc_target_likelihood -1"
    end if
end function

! ------------------------------------------------------------------------------

! ******************************************************************************
! METROPOLIS_HASTINGS TEST METHODS
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
    type(mcmc_sampler) :: mcmc

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
end module
! SPDX-FileCopyrightText: 2022-2026 Jason Christopherson
! SPDX-License-Identifier: MIT
!
! MIT License
!
! Copyright (c) 2022-2026 Jason Christopherson
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module fstats_mcmc_tests
    use iso_fortran_env
    use fstats
    use fortran_test_helper
    use, intrinsic :: ieee_arithmetic
    implicit none

    type, extends(mcmc_target) :: test_mcmc_target
    contains
        procedure, public :: model => tmt_eval
    end type

    type, extends(mcmc_proposal) :: invalid_mcmc_proposal
    contains
        procedure, public :: generate_sample => invalid_proposal
    end type

contains
subroutine invalid_proposal(this, tgt, xc, xp, vc, vp)
    !! Supplies a deterministic NaN log-variance proposal to test rejection.
    class(invalid_mcmc_proposal), intent(inout) :: this
        !! Test proposal object.
    class(mcmc_target), intent(inout) :: tgt
        !! Target, unused by this deterministic proposal.
    real(real64), intent(in), dimension(:) :: xc
        !! Current model parameters.
    real(real64), intent(out), dimension(:) :: xp
        !! Proposed parameters, unchanged from xc.
    real(real64), intent(in) :: vc
        !! Current log variance, unused.
    real(real64), intent(out) :: vp
        !! NaN proposed log variance.

    xp = xc
    vp = ieee_value(0.0d0, ieee_quiet_nan)
end subroutine

function test_mcmc_numerical_robustness() result(rst)
    !! Checks log-variance extremes, the transformed prior, and invalid proposals.
    logical :: rst
        !! True when all deterministic MCMC numerical contracts hold.
    real(real64), dimension(2) :: params, xdata, ydata
    real(real64), allocatable, dimension(:,:) :: chain
    real(real64) :: actual, expected, log_variance, constant
    type(test_mcmc_target) :: target
    type(normal_distribution) :: prior
    type(invalid_mcmc_proposal) :: proposal
    type(mcmc_sampler) :: sampler

    rst = .true.
    call prior%standardize()
    call target%add_parameter(prior)
    call target%add_parameter(prior)
    params = 0.0d0
    xdata = 0.0d0
    ydata = 0.0d0
    constant = log(4.0d0 * acos(0.0d0))
    log_variance = 1000.0d0
    actual = target%likelihood_log_variance(xdata, ydata, params, log_variance)
    expected = -constant - log_variance
    rst = rst .and. abs(actual - expected) < 1.0d-12
    actual = target%likelihood_log_variance(xdata, ydata, params, -log_variance)
    expected = -constant + log_variance
    rst = rst .and. abs(actual - expected) < 1.0d-12
    log_variance = 2.0d0
    actual = target%log_posterior(xdata(:0), ydata(:0), params, log_variance)
    expected = target%evaluate_prior(params) + target%evaluate_variance_prior(exp(log_variance)) + log_variance
    rst = rst .and. abs(actual - expected) < 1.0d-13
    ydata = huge(1.0d0)
    params(2) = -huge(1.0d0)
    log_variance = 2.0d0 * log(huge(1.0d0))
    actual = target%likelihood_log_variance(xdata, ydata, params, log_variance)
    expected = -4.0d0 - constant - log_variance
    rst = rst .and. ieee_is_finite(actual)
    rst = rst .and. abs(actual - expected) < 1.0d-10
    params = 0.0d0
    ydata = 1.0d0
    actual = target%likelihood_log_variance(xdata, ydata, params, -1000.0d0)
    rst = rst .and. .not.ieee_is_finite(actual) .and. actual < 0.0d0
    params = 40.0d0
    rst = rst .and. ieee_is_finite(target%evaluate_prior(params))
    ydata = 0.0d0
    call sampler%sample(xdata, ydata, proposal, target, 3_int32)
    rst = rst .and. sampler%get_accepted_count() == 0
    chain = sampler%get_chain()
    rst = rst .and. all(chain == 0.0d0)
    if (.not.rst) print '(A)', 'TEST FAILED: MCMC numerical robustness'
end function

! ******************************************************************************
! TEST_MCMC_TARGET
! ------------------------------------------------------------------------------
pure subroutine tmt_eval(this, xdata, xc, y)
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
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Local Variables
    integer(int32) :: i
    real(real64) :: v, sm, l, xdata(ndata), ydata(ndata), ymod(ndata), params(2)
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

    ! Evaluate the Gaussian log-density directly to avoid underflow from the
    ! raw PDF in the far tail of the distribution.  The target likelihood is
    ! computed in log-space as a sum of log densities.
    dist%standard_deviation = sqrt(var)
    sm = 0.0d0
    do i = 1, ndata
        dist%mean_value = ymod(i)
        v = -0.5d0 * ((ydata(i) - dist%mean_value) / dist%standard_deviation)**2 &
            - log(dist%standard_deviation * sqrt(2.0d0 * pi))
        sm = sm + v
    end do

    ! Compute the likelihood via the target type
    l = target%likelihood(xdata, ydata, params, var)

    ! Test
    if (.not.assert(sm, l, tol)) then
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
function test_mcmc_parallel_likelihood() result(rst)
    logical :: rst
    integer(int32), parameter :: ndata = 20001
    integer(int32) :: i, mode, variance_index
    real(real64) :: params(2), variances(3), effective_variance, expected, actual
    real(real64) :: xdata(ndata), ydata(ndata), ymod(ndata)
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    type(normal_distribution) :: prior
    type(test_mcmc_target) :: target

    rst = .true.
    params = [2.0d0, 1.0d0]
    variances = [0.25d0, 0.0d0, -1.0d0]
    call target%add_parameter(prior)
    call target%add_parameter(prior)
    do i = 1, ndata
        xdata(i) = real(i - 1, real64) / real(ndata - 1, real64)
    end do
    call target%model(xdata, params, ymod)
    ydata = ymod + 0.125d0 * sin(xdata)

    do mode = 1, 3
        select case (mode)
        case (1)
            target%likelihood_parallel_threshold = 10000
        case (2)
            target%likelihood_parallel_threshold = huge(1_int32)
        case (3)
            target%likelihood_parallel_threshold = 0
        end select
        do variance_index = 1, size(variances)
            effective_variance = variances(variance_index)
            if (effective_variance <= 0.0d0) then
                actual = target%likelihood(xdata, ydata, params, effective_variance)
                if (.not.ieee_is_nan(actual)) rst = .false.
                cycle
            end if
            expected = sum(-(ydata - ymod)**2 / (2.0d0 * effective_variance) &
                - log(sqrt(2.0d0 * pi * effective_variance)))
            actual = target%likelihood(xdata, ydata, params, variances(variance_index))
            if (.not.ieee_is_finite(actual) .or. &
                abs(actual - expected) > 1.0d-10 * max(1.0d0, abs(expected))) rst = .false.
        end do
    end do

    actual = target%likelihood(xdata(:0), ydata(:0), params, 1.0d0)
    if (actual /= 0.0d0) rst = .false.
    actual = target%likelihood(xdata(:21), ydata(:21), params, 0.25d0)
    expected = sum(-(ydata(:21) - ymod(:21))**2 / 0.5d0 - log(sqrt(0.5d0 * pi)))
    if (.not.assert(actual, expected, 1.0d-10)) rst = .false.

    ydata = 1.0d153
    expected = sum(-(ydata - ymod)**2 / 2.0d300 - log(sqrt(2.0d300 * pi)))
    do mode = 1, 2
        target%likelihood_parallel_threshold = 0
        if (mode == 2) target%likelihood_parallel_threshold = huge(1_int32)
        actual = target%likelihood(xdata, ydata, params, 1.0d300)
        if (.not.ieee_is_finite(actual) .or. &
            abs(actual - expected) > 1.0d-10 * abs(expected)) rst = .false.
    end do
    if (.not.rst) print '(A)', 'TEST FAILED: test_mcmc_parallel_likelihood'
end function

! ------------------------------------------------------------------------------
function test_mcmc_sample_chains() result(rst)
    logical :: rst
    integer(int32), parameter :: nchains = 4, niter = 205, ndata = 21
    integer(int32) :: chain_index, i
    real(real64) :: xdata(ndata), ydata(ndata)
    real(real64), allocatable :: chain(:,:)
    type(normal_distribution) :: prior
    type(test_mcmc_target) :: targets(nchains)
    type(mcmc_proposal) :: proposals(nchains)
    type(mcmc_sampler) :: samplers(nchains)

    rst = .true.
    do i = 1, ndata
        xdata(i) = real(i - 1, real64) / real(ndata - 1, real64)
    end do
    ydata = 2.0d0 * xdata + 1.0d0
    prior%standard_deviation = 1.0d0
    do chain_index = 1, nchains
        prior%mean_value = real(chain_index, real64)
        call targets(chain_index)%add_parameter(prior)
        prior%mean_value = -real(chain_index, real64)
        call targets(chain_index)%add_parameter(prior)
        targets(chain_index)%likelihood_parallel_threshold = 0
        call proposals(chain_index)%set_scale(0.01d0 * real(chain_index, real64))
    end do

    call sample_chains(samplers, xdata, ydata, proposals, targets, niter)
    do chain_index = 1, nchains
        if (samplers(chain_index)%get_chain_length() /= niter) rst = .false.
        if (samplers(chain_index)%get_state_variable_count() /= 3) rst = .false.
        if (samplers(chain_index)%get_accepted_count() < 0 .or. &
            samplers(chain_index)%get_accepted_count() >= niter) rst = .false.
        chain = samplers(chain_index)%get_chain()
        if (.not.all(ieee_is_finite(chain))) rst = .false.
        if (.not.assert(chain(1,:), &
            [real(chain_index, real64), -real(chain_index, real64), 0.0d0])) rst = .false.
    end do

    call sample_chains(samplers, xdata, ydata, proposals, targets, 1_int32)
    do chain_index = 1, nchains
        if (samplers(chain_index)%get_chain_length() /= niter + 1) rst = .false.
        if (samplers(chain_index)%get_accepted_count() /= 0) rst = .false.
        chain = samplers(chain_index)%get_chain()
        if (.not.assert(chain(1,:), chain(niter + 1,:))) rst = .false.
    end do

    targets(1)%likelihood_parallel_threshold = 10000
    call sample_chains(samplers(:1), xdata, ydata, proposals(:1), targets(:1))
    if (samplers(1)%get_chain_length() /= niter + 1 + 10000) rst = .false.
    call sample_chains(samplers(:0), xdata, ydata, proposals(:0), targets(:0))
    if (.not.rst) print '(A)', 'TEST FAILED: test_mcmc_sample_chains'
end function

! ------------------------------------------------------------------------------
end module
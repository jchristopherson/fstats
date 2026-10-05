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

module fstats_bootstrap_tests
    use iso_fortran_env
    use fstats
    use fortran_test_helper
    implicit none
contains
    function test_bootstrap_1() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: tol = 5.0d-2
        real(real64), parameter :: ci_upper = 0.746d0
        real(real64), parameter :: ci_lower = 0.583d0
        real(real64), parameter :: std_err = 0.042d0
        integer(int32), parameter :: nsamples = 1000

        ! NOTES:
        ! - The loose tolerance accounts for the precision with which
        ! the answers were computed and to account for the different
        ! resampling distributions that may arrise from different calls
        ! to the bootstrap routine itself.
        !
        ! - Solutions computed via JASP
        !
        ! - Original data set was randomly generated

        ! Variables
        integer(int32), parameter :: n = 30
        real(real64) :: x(n), mu
        procedure(bootstrap_statistic_routine), pointer :: fcn
        type(bootstrap_statistics) :: z

        ! Initialization
        rst = .true.
        x = [0.986840296d0, 0.932169009d0, 0.870384823d0, 0.532368509d0, &
            0.707833926d0, 0.447913731d0, 0.390929201d0, 0.770973409d0, &
            0.958743382d0, 0.957561919d0, 0.470107632d0, 0.768371184d0, &
            0.725434044d0, 0.311942491d0, 0.673420327d0, 0.897531681d0, &
            0.669168983d0, 0.544250163d0, 0.211680387d0, 0.768569119d0, &
            0.988756414d0, 0.214287163d0, 0.579534890d0, 0.770211065d0, &
            0.307955851d0, 0.509150720d0, 0.665120628d0, 0.817904438d0, &
            0.896080063d0, 0.599850776d0]
        fcn => mean
        mu = mean(x)

        ! Test
        z = bootstrap(fcn, x, nsamples = nsamples)

        if (.not.assert(z%statistic_value, mu, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_bootstrap_1 -1"
        end if

        if (.not.assert(z%upper_confidence_interval, ci_upper, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_bootstrap_1 -2"
        end if

        if (.not.assert(z%lower_confidence_interval, ci_lower, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_bootstrap_1 -3"
        end if

        if (.not.assert(z%standard_error, std_err, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_bootstrap_1 -4"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_random_resample_with_replacement() result(rst)
        logical :: rst

        real(real64) :: x(5), xn(5)

        rst = .true.
        x = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0]
        call random_resample_with_replacement(x, xn)

        if (size(xn) /= size(x)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_random_resample_with_replacement -1"
        end if

        if (.not.all(xn >= minval(x) .and. xn <= maxval(x))) then
            rst = .false.
            print "(A)", "TEST FAILED: test_random_resample_with_replacement -2"
        end if

        if (.not.any(xn == x(1)) .and. .not.any(xn == x(2)) .and. &
            .not.any(xn == x(3)) .and. .not.any(xn == x(4)) .and. &
            .not.any(xn == x(5))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_random_resample_with_replacement -3"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_bootstrap_2() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: tol = 5.0d-2
        real(real64), parameter :: ci_upper = 0.746d0
        real(real64), parameter :: ci_lower = 0.583d0
        real(real64), parameter :: std_err = 0.042d0
        integer(int32), parameter :: nsamples = 1000

        ! NOTES:
        ! - The loose tolerance accounts for the precision with which
        ! the answers were computed and to account for the different
        ! resampling distributions that may arrise from different calls
        ! to the bootstrap routine itself.
        !
        ! - Solutions computed via JASP
        !
        ! - Original data set was randomly generated

        ! Variables
        integer(int32), parameter :: n = 30
        real(real64) :: x(n), mu
        procedure(bootstrap_statistic_routine), pointer :: fcn
        procedure(bootstrap_resampling_routine), pointer :: sampler
        type(bootstrap_statistics) :: z

        ! Initialization
        rst = .true.
        x = [0.986840296d0, 0.932169009d0, 0.870384823d0, 0.532368509d0, &
            0.707833926d0, 0.447913731d0, 0.390929201d0, 0.770973409d0, &
            0.958743382d0, 0.957561919d0, 0.470107632d0, 0.768371184d0, &
            0.725434044d0, 0.311942491d0, 0.673420327d0, 0.897531681d0, &
            0.669168983d0, 0.544250163d0, 0.211680387d0, 0.768569119d0, &
            0.988756414d0, 0.214287163d0, 0.579534890d0, 0.770211065d0, &
            0.307955851d0, 0.509150720d0, 0.665120628d0, 0.817904438d0, &
            0.896080063d0, 0.599850776d0]
        fcn => mean
        sampler => random_resample_with_replacement
        mu = mean(x)

        ! Test
        z = bootstrap(fcn, x, method = sampler, nsamples = nsamples)

        if (.not.assert(z%statistic_value, mu, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_bootstrap_2 -1"
        end if

        if (.not.assert(z%upper_confidence_interval, ci_upper, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_bootstrap_2 -2"
        end if

        if (.not.assert(z%lower_confidence_interval, ci_lower, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_bootstrap_2 -3"
        end if

        if (.not.assert(z%standard_error, std_err, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_bootstrap_2 -4"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_bootstrap_small_nsamples() result(rst)
        logical :: rst
        real(real64) :: x(3)
        procedure(bootstrap_statistic_routine), pointer :: fcn
        procedure(bootstrap_resampling_routine), pointer :: sampler
        type(bootstrap_statistics) :: z

        x = [1.0d0, 2.0d0, 3.0d0]
        fcn => mean
        sampler => shift_resample
        z = bootstrap(fcn, x, method = sampler, nsamples = 3)

        rst = z%statistic_value == 2.0d0 .and. &
            all(z%population == 12.0d0) .and. &
            z%lower_confidence_interval == 12.0d0 .and. &
            z%upper_confidence_interval == 12.0d0
        if (.not.rst) print "(A)", "TEST FAILED: test_bootstrap_small_nsamples"
    end function

! ------------------------------------------------------------------------------
    function test_structure_preserving_resamplers() result(rst)
        logical :: rst
        real(real64) :: paired(8), paired_sample(8)
        real(real64) :: clustered(8), cluster_sample(8)
        real(real64) :: series(12), block_sample(12)
        integer(int32) :: i, block_size

        rst = .true.
        paired = [1.0d0, 101.0d0, 2.0d0, 102.0d0, &
            3.0d0, 103.0d0, 4.0d0, 104.0d0]
        call random_resample_paired(paired, paired_sample)
        do i = 1, size(paired_sample), 2
            if (paired_sample(i + 1) - paired_sample(i) /= 100.0d0) rst = .false.
        end do

        clustered = [1.0d0, 10.0d0, 1.0d0, 11.0d0, &
            2.0d0, 20.0d0, 2.0d0, 21.0d0]
        call random_resample_clusters(clustered, cluster_sample)
        if (cluster_sample(1) /= cluster_sample(3)) rst = .false.
        if (cluster_sample(5) /= cluster_sample(7)) rst = .false.
        do i = 1, size(cluster_sample), 2
            if (cluster_sample(i) == 1.0d0) then
                if (cluster_sample(i + 1) < 10.0d0 .or. &
                    cluster_sample(i + 1) > 11.0d0) rst = .false.
            else if (cluster_sample(i) == 2.0d0) then
                if (cluster_sample(i + 1) < 20.0d0 .or. &
                    cluster_sample(i + 1) > 21.0d0) rst = .false.
            else
                rst = .false.
            end if
        end do

        series = [(real(i, real64), i = 1, size(series))]
        call random_resample_circular_blocks(series, block_sample)
        block_size = max(2_int32, int(sqrt(real(size(series), real64)), int32))
        if (any(block_sample < 1.0d0 .or. block_sample > real(size(series), real64))) &
            rst = .false.
        do i = 2, size(block_sample)
            if (mod(i - 1, block_size) /= 0) then
                if (block_sample(i) /= &
                    real(modulo(int(block_sample(i - 1), int32), &
                    size(series)) + 1, real64)) rst = .false.
            end if
        end do

        if (.not.rst) print "(A)", "TEST FAILED: test_structure_preserving_resamplers"
    end function

! ------------------------------------------------------------------------------
    subroutine shift_resample(x, xn)
        real(real64), intent(in) :: x(:)
        real(real64), intent(out) :: xn(size(x))

        xn = x + 10.0d0
    end subroutine

! ------------------------------------------------------------------------------
end module
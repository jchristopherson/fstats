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
        real(real64), parameter :: tol = 1.0d-2
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

        real(real64), parameter :: tol = 1.0d-12
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
        real(real64), parameter :: tol = 1.0d-2
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
        z = bootstrap(fcn, x, nsamples = nsamples)

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
end module
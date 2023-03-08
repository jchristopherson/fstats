module fstats_distribution_tests
    use iso_fortran_env
    use fstats
    use fstats_test_helper
    implicit none
contains
! ------------------------------------------------------------------------------
    function t_distribution_test_1() result(rst)
        ! Arguments
        logical :: rst

        ! Variables
        integer(int32), parameter :: n = 9
        real(real64) :: x(n), pdf(n), cdf(n), pdfans(n), cdfans(n)
        type(t_distribution) :: dist

        ! Initialization - solution computed via Excel
        rst = .true.
        dist%dof = 20.0d0
        x = [-1.0d0, -0.75d0, -0.5d0, -0.25d0, 0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0]
        pdfans = [0.23604564912670d0, 0.29444316943118d0, 0.34580861238374d0, &
            0.38129013769196d0, 0.3939885857114d0, 0.3812901376920d0, &
            0.3458086123837d0, 0.2944431694312d0, 0.2360456491267d0]
        cdfans = [0.16462828858586d0, 0.23099351240633d0, 0.31126592114051d0, &
            0.40256865848010d0, 0.5d0, 0.5974313415199d0, 0.6887340788595d0, &
            0.7690064875937d0, 0.8353717114141d0]

        ! Tests
        pdf = dist%pdf(x)
        if (.not.is_equal(pdf, pdfans)) then
            rst = .false.
            print '(A)', "TEST FAILED: T-Distribution PDF test."
        end if

        cdf = dist%cdf(x)
        if (.not.is_equal(cdf, cdfans)) then
            rst = .false.
            print '(A)', "TEST FAILED: T-Distribution CDF test."
        end if
    end function

! ------------------------------------------------------------------------------
    function normal_distribution_test_1() result(rst)
        ! Arguments
        logical :: rst

        ! Variables
        integer(int32), parameter :: n = 9
        real(real64) :: x(n), pdf(n), cdf(n), pdfans(n), cdfans(n)
        type(normal_distribution) :: dist

        ! Initialization
        rst = .true.
        call dist%standardize()
        x = [-1.0d0, -0.75d0, -0.5d0, -0.25d0, 0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0]
        pdfans = [0.24197072451914d0, 0.30113743215480d0, 0.35206532676430d0, &
            0.38666811680285d0, 0.39894228040143d0, 0.38666811680285d0, &
            0.35206532676430d0, 0.30113743215480d0, 0.24197072451914d0]
        cdfans = [0.15865525393146d0, 0.22662735237687d0, 0.30853753872599d0, &
            0.40129367431708d0, 0.5d0, 0.59870632568292d0, 0.69146246127401d0, &
            0.77337264762313d0, 0.84134474606854d0]

        ! Tests
        pdf = dist%pdf(x)
        if (.not.is_equal(pdf, pdfans)) then
            rst = .false.
            print '(A)', "TEST FAILED: Normal Distribution PDF test."
        end if

        cdf = dist%cdf(x)
        if (.not.is_equal(cdf, cdfans)) then
            rst = .false.
            print '(A)', "TEST FAILED: Normal Distribution CDF test."
        end if
    end function

! ------------------------------------------------------------------------------
    function f_distribution_test_1() result(rst)
        ! Arguments
        logical :: rst

        ! Variables
        integer(int32), parameter :: n = 9
        real(real64) :: x(n), pdf(n), cdf(n), pdfans(n), cdfans(n)
        type(f_distribution) :: dist

        ! Initialization
        rst = .true.
        dist%d1 = 8.0d0
        dist%d2 = 10.0d0
        x = [0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0, 1.25d0, 1.5d0, 1.75d0, 2.0d0]
        pdfans = [0.0d0, 0.347301605446d0, 0.693866102230d0, 0.704079866409d0, &
            0.578183153344d0, 0.437500000000d0, 0.320617799490d0, &
            0.232664942711d0, 0.168984874907d0]
        cdfans = [0.0d0, 0.030656411942d0, 0.168986926001d0, 0.348632991314d0, &
            0.510052693677d0, 0.636718750000d0, 0.730868976686d0, &
            0.799460940309d0, 0.849226439763d0]

        ! Tests
        pdf = dist%pdf(x)
        if (.not.is_equal(pdf, pdfans)) then
            rst = .false.
            print '(A)', "TEST FAILED: F Distribution PDF test."
        end if

        cdf = dist%cdf(x)
        if (.not.is_equal(cdf, cdfans)) then
            rst = .false.
            print '(A)', "TEST FAILED: F Distribution CDF test."
        end if
    end function

! ------------------------------------------------------------------------------
    function chi_squared_distribution_test_1() result(rst)
        ! Arguments
        logical :: rst

        ! Variables
        integer(int32), parameter :: dof = 5
        integer(int32), parameter :: n = 9
        real(real64) :: x(n), pdf(n), cdf(n), pdfans(n), cdfans(n)
        type(chi_squared_distribution) :: dist

        ! Initialization
        rst = .true.
        dist%dof = dof
        x = [0.0d0, 0.5d0, 1.0d0, 1.5d0, 2.0d0, 2.5d0, 3.0d0, 3.5d0, 4.0d0]
        pdfans = [0.0d0, 0.036615940789d0, 0.080656908173d0, 0.115399742104d0, &
            0.138369165807d0, 0.150601993890d0, 0.154180329804d0, &
            0.151312753472d0, 0.143975910702d0]
        cdfans = [0.0d0, 0.007876706767d0, 0.037434226753d0, 0.086930185456d0, &
            0.150854963915d0, 0.223504928877d0, 0.300014164121d0, &
            0.376612372250d0, 0.450584048647d0]

        ! Tests
        pdf = dist%pdf(x)
        if (.not.is_equal(pdf, pdfans)) then
            rst = .false.
            print '(A)', "TEST FAILED: Chi-squared distribution PDF test."
        end if

        cdf = dist%cdf(x)
        if (.not.is_equal(cdf, cdfans)) then
            rst = .false.
            print '(A)', "TEST FAILED: Chi-squared distribution CDF test."
        end if
    end function

! ------------------------------------------------------------------------------
end module
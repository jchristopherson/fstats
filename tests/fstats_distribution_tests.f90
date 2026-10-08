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

module fstats_distribution_tests
    use iso_fortran_env
    use fstats
    use fstats_test_helper
    use ieee_arithmetic
    implicit none
contains
    function test_distribution_numerical_robustness() result(rst)
        !! Checks extreme tails, log densities, discrete boundaries, and invalid inputs.
        logical :: rst
            !! True when every numerical stress case passes.
        type(normal_distribution) :: normal
        type(log_normal_distribution) :: lognormal
        type(t_distribution) :: student
        type(binomial_distribution) :: binomial
        type(poisson_distribution) :: poisson
        type(chi_squared_distribution) :: chisq
        type(f_distribution) :: fisher
        real(real64) :: value, nan

        rst = .true.
        nan = ieee_value(0.0d0, ieee_quiet_nan)
        call normal%standardize()
        rst = rst .and. normal%pdf(40.0d0) == 0.0d0
        rst = rst .and. ieee_is_finite(normal%log_pdf(40.0d0))
        value = 7.619853024160526d-24
        rst = rst .and. abs(normal%survival(10.0d0) / value - 1.0d0) < 1.0d-13
        rst = rst .and. abs(normal%cdf(-10.0d0) / value - 1.0d0) < 1.0d-13
        rst = rst .and. abs(normal%log_survival(40.0d0) + 804.6084420137538d0) < 1.0d-11
        rst = rst .and. normal%log_cdf(-40.0d0) == normal%log_survival(40.0d0)
        normal%mean_value = -1.0d308
        normal%standard_deviation = 1.0d308
        rst = rst .and. abs(normal%log_pdf(1.0d308) + 2.0d0 + log(1.0d308) &
            + 0.5d0 * log(4.0d0 * acos(0.0d0))) < 1.0d-12
        call normal%standardize()
        lognormal%mean_value = 2.0d0
        lognormal%standard_deviation = 0.5d0
        rst = rst .and. abs(lognormal%cdf(exp(2.0d0)) - 0.5d0) < 1.0d-14
        rst = rst .and. lognormal%pdf(0.0d0) == 0.0d0
        rst = rst .and. lognormal%survival(-1.0d0) == 1.0d0
        rst = rst .and. ieee_is_nan(lognormal%log_pdf(nan))
        lognormal%mean_value = 0.0d0
        lognormal%standard_deviation = 1.0d-10
        rst = rst .and. abs(lognormal%variance() / 1.0d-20 - 1.0d0) < 1.0d-14
        student%dof = 1000.0d0
        rst = rst .and. abs(student%pdf(0.0d0) - 0.3988425573138582d0) < 1.0d-12
        rst = rst .and. student%cdf(-20.0d0) > 0.0d0
        binomial%n = 1000
        binomial%p = 0.5d0
        rst = rst .and. abs(binomial%pdf(500.0d0) - 0.02522501817836080d0) < 1.0d-12
        rst = rst .and. binomial%pdf(500.5d0) == 0.0d0
        rst = rst .and. binomial%cdf(1000.0d0) == 1.0d0
        binomial%p = 0.0d0
        rst = rst .and. binomial%pdf(0.0d0) == 1.0d0
        poisson%occrence_rate = 1000.0d0
        rst = rst .and. abs(poisson%pdf(1000.0d0) - 0.01261461134872150d0) < 1.0d-12
        rst = rst .and. abs(poisson%cdf(1000.0d0) + poisson%survival(1000.0d0) - 1.0d0) < 1.0d-13
        poisson%occrence_rate = 0.0d0
        rst = rst .and. poisson%pdf(0.0d0) == 1.0d0
        chisq%dof = 2
        rst = rst .and. abs(chisq%survival(100.0d0) / exp(-50.0d0) - 1.0d0) < 1.0d-13
        fisher%d1 = 2.0d0
        fisher%d2 = 2.0d0
        rst = rst .and. abs(fisher%survival(1.0d20) / 1.0d-20 - 1.0d0) < 1.0d-13
        normal%standard_deviation = 0.0d0
        rst = rst .and. ieee_is_nan(normal%pdf(0.0d0))
        chisq%dof = 2
        rst = rst .and. chisq%survival(2000.0d0) == 0.0d0
        rst = rst .and. abs(chisq%log_survival(2000.0d0) + 1000.0d0) < 1.0d-12
        poisson%occrence_rate = 1000.0d0
        rst = rst .and. abs(poisson%log_cdf(0.0d0) + 1000.0d0) < 1.0d-12
        if (.not.rst) print '(A)', 'TEST FAILED: distribution numerical robustness'
    end function

    function test_multivariate_normal_numerical_robustness() result(rst)
        !! Checks Cholesky log densities with overflowing/underflowing determinants.
        logical :: rst
            !! True when both covariance scales give the analytic density.
        type(multivariate_normal_distribution) :: normal
        real(real64), dimension(2,2) :: covariance_values
        real(real64), dimension(2) :: point
        real(real64) :: expected

        rst = .true.
        covariance_values = 0.0d0
        covariance_values(1,1) = 1.0d200
        covariance_values(2,2) = 1.0d200
        point = 1.0d100
        call normal%initialize([0.0d0, 0.0d0], covariance_values)
        expected = -1.0d0 - log(4.0d0 * acos(0.0d0)) - log(1.0d200)
        rst = rst .and. abs(normal%log_pdf(point) - expected) < 1.0d-12
        covariance_values(1,1) = 1.0d-200
        covariance_values(2,2) = 1.0d-200
        point = 1.0d-100
        call normal%initialize([0.0d0, 0.0d0], covariance_values)
        expected = -1.0d0 - log(4.0d0 * acos(0.0d0)) - log(1.0d-200)
        rst = rst .and. abs(normal%log_pdf(point) - expected) < 1.0d-12
        if (.not.rst) print '(A)', 'TEST FAILED: multivariate-normal numerical robustness'
    end function

    function test_special_function_numerical_robustness() result(rst)
        !! Checks log tails, small-argument helpers, invalid domains, and iteration limits.
        logical :: rst
            !! True when all special-function numerical contracts hold.

        rst = abs(log_regularized_gamma_upper(1.0d0, 1000.0d0) + 1000.0d0) < 1.0d-12
        rst = rst .and. abs(log_regularized_beta(1000.0d0, 1.0d0, 0.1d0) &
            - 1000.0d0 * log(0.1d0)) < 1.0d-10
        rst = rst .and. abs(stable_log1p(1.0d-20) / 1.0d-20 - 1.0d0) < 1.0d-14
        rst = rst .and. abs(stable_expm1(1.0d-20) / 1.0d-20 - 1.0d0) < 1.0d-14
        rst = rst .and. ieee_is_nan(regularized_beta(2.0d0, 3.0d0, 0.4d0, 0_int32))
        rst = rst .and. ieee_is_nan(regularized_gamma_lower(2.0d0, 1.0d0, 0_int32))
        rst = rst .and. ieee_is_nan(regularized_gamma_upper(2.0d0, 4.0d0, 0_int32))
        rst = rst .and. ieee_is_nan(regularized_beta(-1.0d0, 1.0d0, 0.5d0))
        if (.not.rst) print '(A)', 'TEST FAILED: special-function numerical robustness'
    end function

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
    function binomial_distribution_test_1() result(rst)
        ! Arguments
        logical :: rst

        ! Variables
        integer(int32), parameter :: n = 100
        real(real64), parameter :: p = 0.25d0
        integer(int32), parameter :: npts = 9
        real(real64) :: x(npts), pmf(npts), cdf(npts), pmfans(npts), cdfans(npts)
        type(binomial_distribution) :: dist

        ! Initialization
        rst = .true.
        dist%n = n
        dist%p = p
        x = [1.0d1, 2.0d1, 3.0d1, 4.0d1, 5.0d1, 6.0d1, 7.0d1, 8.0d1, 9.0d1]
        pmfans = [ &
            9.40196486279606d-05, 4.93006403376772d-02, 4.57538076104868d-02, &
            3.62626791645609d-04, 4.50731087508638d-08, 1.04000348155055d-13, &
            3.76337116402577d-21, 1.16299347183680d-30, 6.36089528619433d-43]
        cdfans = [ &
            0.000137100563168d0, 0.148831050442992d0, 0.896212761043913d0, &
            0.999676034583683d0, 0.999999978688084d0, 0.999999999999971d0, &
            1.0d0, 1.0d0, 1.0d0]

        ! Tests
        pmf = dist%pdf(x)
        if (.not.is_equal(pmf, pmfans)) then
            rst = .false.
            print '(A)', "TEST FAILED: Binomial distribution PMF test."
        end if

        cdf = dist%cdf(x)
        if (.not.is_equal(cdf, cdfans)) then
            rst = .false.
            print '(A)', "TEST FAILED: Binomial distribution CDF test."
        end if
    end function

! ------------------------------------------------------------------------------
    function test_standardized_variable() result(rst)
        ! Arguments
        logical :: rst

        ! Variables
        real(real64), parameter :: tol = 1.0d-2
        real(real64), parameter :: alpha1 = 0.975d0
        real(real64), parameter :: alpha2 = 0.2d0
        real(real64), parameter :: ans1 = 1.96d0
        real(real64), parameter :: ans2 = 0.84d0
        real(real64) :: z1, z2
        type(normal_distribution) :: dist

        ! Initialization
        rst = .true.
        call dist%standardize()

        ! Test 1
        z1 = dist%standardized_variable(alpha1)
        if (.not.is_equal(z1, ans1, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: Standardized variable test -1."
        end if

        ! Test 2
        z2 = dist%standardized_variable(0.8d0)
        if (.not.is_equal(z2, ans2, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: Standardized variable test -2."
        end if
    end function

! ------------------------------------------------------------------------------
    function test_multivariate_normal_distribution() result(rst)
        use linalg, only : mtx_inverse, det
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64) :: x(2), mu(2), rho, s1, s2, sigma(2, 2), arg, ans, phi, &
            dsig, inv(2, 2)
        type(multivariate_normal_distribution) :: dist

        ! Initialization
        rst = .true.
        call random_number(x)
        call random_number(mu)
        call random_number(rho)
        call random_number(s1)
        call random_number(s2)
        sigma = reshape([s2**2, -rho * s1 * s2, -rho * s1 * s2, s1**2], [2, 2])
        call dist%initialize(mu, sigma)

        ! Compute the actual solution
        inv = mtx_inverse(sigma)
        arg = -0.5d0 * dot_product(x - mu, matmul(inv, x - mu))
        dsig = det(sigma)
        ans = exp(arg) / sqrt((2.0d0 * pi)**2 * dsig)

        ! Test
        phi = dist%pdf(x)
        if (.not.is_equal(phi, ans, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_multivariate_normal_distribution -1"
        end if
    end function

! ------------------------------------------------------------------------------
end module
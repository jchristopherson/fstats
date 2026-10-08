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

module fstats_distributions
    use iso_fortran_env
    use ieee_arithmetic
    use fstats_special_functions
    use fstats_helper_routines
    use fstats_errors
    implicit none
    private
    public :: distribution
    public :: distribution_function
    public :: distribution_property
    public :: distribution_recenter
    public :: t_distribution
    public :: normal_distribution
    public :: f_distribution
    public :: chi_squared_distribution
    public :: binomial_distribution
    public :: multivariate_distribution
    public :: multivariate_distribution_function
    public :: multivariate_normal_distribution
    public :: log_normal_distribution
    public :: poisson_distribution

    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    type, abstract :: distribution
        !! Defines a probability distribution.
    contains
        procedure(distribution_function), public, deferred, pass :: pdf
        procedure(distribution_function), public, deferred, pass :: cdf
        procedure(distribution_property), public, deferred, pass :: mean
        procedure(distribution_property), public, deferred, pass :: median
        procedure(distribution_property), public, deferred, pass :: mode
        procedure(distribution_property), public, deferred, pass :: variance
        procedure(distribution_recenter), public, deferred, pass :: recenter
        procedure, public :: standardized_variable => dist_std_var
        procedure, public :: defined_range => dist_defined_range
        procedure, public :: log_pdf => dist_log_pdf
            !! Natural log density (or log probability mass for discrete laws).
        procedure, public :: survival => dist_survival
            !! Upper-tail probability P(X > x), evaluated directly for built-ins.
        procedure, public :: log_cdf => dist_log_cdf
            !! Natural log of P(X <= x).
        procedure, public :: log_survival => dist_log_survival
            !! Natural log of P(X > x).
    end type

    interface
        pure elemental function distribution_function(this, x) result(rst)
            !! Defines the interface for a probability distribution function.
            use iso_fortran_env, only : real64
            import distribution
            class(distribution), intent(in) :: this
                !! The distribution object.
            real(real64), intent(in) :: x
                !! The value at which to evaluate the function.
            real(real64) :: rst
                !! The value of the function.
        end function

        pure function distribution_property(this) result(rst)
            !! Computes the value of a distribution property.
            use iso_fortran_env, only : real64
            import distribution
            class(distribution), intent(in) :: this
                !! The distribution object.
            real(real64) :: rst
                !! The property value.
        end function

        subroutine distribution_recenter(this, x)
            !! Recenters the distribution about the supplied value.
            use iso_fortran_env, only : real64
            import distribution
            class(distribution), intent(inout) :: this
                !! The distribution object.
            real(real64), intent(in) :: x
                !! The value about which to recenter.
        end subroutine
    end interface

! ------------------------------------------------------------------------------
    type, extends(distribution) :: t_distribution
        !! Defines Student's T-Distribution.
        real(real64) :: dof
            !! The number of degrees of freedom.
    contains
        procedure, public :: pdf => td_pdf
        procedure, public :: cdf => td_cdf
        procedure, public :: mean => td_mean
        procedure, public :: median => td_median
        procedure, public :: mode => td_mode
        procedure, public :: variance => td_variance
        procedure, public :: recenter => td_recenter
    end type
! ------------------------------------------------------------------------------
    type, extends(distribution) :: normal_distribution
        !! Defines a normal distribution.
        real(real64) :: standard_deviation
            !! The standard deviation of the distribution.
        real(real64) :: mean_value
            !! The mean value of the distribution.
    contains
        procedure, public :: pdf => nd_pdf
        procedure, public :: cdf => nd_cdf
        procedure, public :: log_pdf => nd_log_pdf
        procedure, public :: survival => nd_survival
        procedure, public :: log_cdf => nd_log_cdf
        procedure, public :: log_survival => nd_log_survival
        procedure, public :: mean => nd_mean
        procedure, public :: median => nd_median
        procedure, public :: mode => nd_mode
        procedure, public :: variance => nd_variance
        procedure, public :: standardize => nd_standardize
        procedure, public :: recenter => nd_recenter
    end type

! ------------------------------------------------------------------------------
    type, extends(distribution) :: f_distribution
        !! Defines an F-distribution.
        real(real64) :: d1
            !! The measure of degrees of freedom for the first data set.
        real(real64) :: d2
            !! The measure of degrees of freedom for the second data set.
    contains
        procedure, public :: pdf => fd_pdf
        procedure, public :: cdf => fd_cdf
        procedure, public :: mean => fd_mean
        procedure, public :: median => fd_median
        procedure, public :: mode => fd_mode
        procedure, public :: variance => fd_variance
        procedure, public :: defined_range => fd_range
        procedure, public :: recenter => fd_recenter
    end type

! ------------------------------------------------------------------------------
    type, extends(distribution) :: chi_squared_distribution
        !! Defines a Chi-squared distribution.
        integer(int32) :: dof
            !! The number of degrees of freedom.
    contains
        procedure, public :: pdf => cs_pdf
        procedure, public :: cdf => cs_cdf
        procedure, public :: mean => cs_mean
        procedure, public :: median => cs_median
        procedure, public :: mode => cs_mode
        procedure, public :: variance => cs_variance
        procedure, public :: defined_range => cs_range
        procedure, public :: recenter => cs_recenter
    end type

! ------------------------------------------------------------------------------
    type, extends(distribution) :: binomial_distribution
        !! Defines a binomial distribution.  The binomial distribution describes
        !! the probability p of getting k successes in n independent trials.
        integer(int32) :: n
            !! The number of independent trials.
        real(real64) :: p
            !! The success probability for each trial.  This parameter must
            !! exist on the set [0, 1].
    contains
        procedure, public :: pdf => bd_pdf
        procedure, public :: cdf => bd_cdf
        procedure, public :: mean => bd_mean
        procedure, public :: median => bd_median
        procedure, public :: mode => bd_mode
        procedure, public :: variance => bd_variance
        procedure, public :: defined_range => bd_range
        procedure, public :: recenter => bd_recenter
    end type

! ------------------------------------------------------------------------------
    type, extends(distribution) :: log_normal_distribution
        !! Defines a normal distribution.
        real(real64) :: standard_deviation
            !! The standard deviation of the distribution.
        real(real64) :: mean_value
            !! The mean value of the distribution.
    contains
        procedure, public :: pdf => lnd_pdf
        procedure, public :: cdf => lnd_cdf
        procedure, public :: log_pdf => lnd_log_pdf
        procedure, public :: survival => lnd_survival
        procedure, public :: log_cdf => lnd_log_cdf
        procedure, public :: log_survival => lnd_log_survival
        procedure, public :: mean => lnd_mean
        procedure, public :: median => lnd_median
        procedure, public :: mode => lnd_mode
        procedure, public :: variance => lnd_variance
        procedure, public :: defined_range => lnd_range
        procedure, public :: recenter => lnd_recenter
    end type

! ------------------------------------------------------------------------------
    type, extends(distribution) :: poisson_distribution
        !! Defines a Poisson distribution.
        real(real64), public :: occrence_rate
            !! The rate of occurrences.
    contains
        procedure, public :: pdf => pd_pdf
        procedure, public :: cdf => pd_cdf
        procedure, public :: mean => pd_mean
        procedure, public :: median => pd_median
        procedure, public :: mode => pd_mode
        procedure, public :: variance => pd_variance
        procedure, public :: recenter => pd_recenter
        procedure, public :: defined_range => pd_range
    end type

! ******************************************************************************
! MULTIVARIATE DISTRIBUTIONS
! ------------------------------------------------------------------------------
    type, abstract :: multivariate_distribution
        !! Defines a multivariate probability distribution.
    contains
        procedure(multivariate_distribution_function), deferred, pass :: pdf
            !! Computes the probability density function.
        procedure, public :: log_pdf => mvd_log_pdf
            !! Computes the log density; custom laws should override the PDF fallback.
    end type

    interface
        pure function multivariate_distribution_function(this, x) result(rst)
            !! Defines an interface for a multivariate probability distribution
            !! function.
            use iso_fortran_env, only : real64
            import multivariate_distribution
            class(multivariate_distribution), intent(in) :: this
                !! The distribution object.
            real(real64), intent(in), dimension(:) :: x
                !! The values at which to evaluate the function.
            real(real64) :: rst
                !! The value of the function.
        end function
    end interface

! ------------------------------------------------------------------------------
    type, extends(multivariate_distribution) :: multivariate_normal_distribution
        !! Defines a multivariate normal (Gaussian) distribution.
        real(real64), private, allocatable, dimension(:) :: m_means
            !! An N-element array of mean values.
        real(real64), private, allocatable, dimension(:,:) :: m_cov
            !! The N-by-N covariance matrix.  This matrix must be 
            !! positive-definite.
        real(real64), private, allocatable, dimension(:,:) :: m_cholesky
            !! The N-by-N Cholesky factored form (lower) of the covariance
            !! matrix.
        real(real64), private :: m_logCovDet
            !! Natural logarithm of the covariance determinant.
    contains
        procedure, public :: initialize => mvnd_init
        procedure, public :: pdf => mvnd_pdf
        procedure, public :: log_pdf => mvnd_log_pdf
        procedure, public :: get_means => mvnd_get_means
        procedure, public :: set_means => mvnd_update_mean
        procedure, public :: get_covariance => mvnd_get_covariance
        procedure, public :: get_cholesky_factored_matrix => mvnd_get_cholesky
    end type

contains
pure elemental function log_one_plus(value) result(rst)
    !! Computes log(1 + value), retaining accuracy for small values.
    real(real64), intent(in) :: value
        !! Argument; values below -1 produce NaN.
    real(real64) :: rst
        !! Logarithm; -1 produces negative infinity.

    rst = stable_log1p(value)
end function

pure elemental function softplus(value) result(rst)
    !! Computes log(1 + exp(value)) without overflowing the exponential.
    real(real64), intent(in) :: value
        !! Exponent argument.
    real(real64) :: rst
        !! Softplus value; NaN inputs propagate.

    rst = max(value, 0.0d0) + log_one_plus(exp(-abs(value)))
end function

pure elemental function probability_log(value) result(rst)
    !! Computes a log density or probability with explicit zero handling.
    real(real64), intent(in) :: value
        !! Nonnegative density or probability.
    real(real64) :: rst
        !! Natural logarithm; zero gives -infinity, negative or NaN gives NaN.

    if (value == 0.0d0) then
        rst = ieee_value(0.0d0, ieee_negative_inf)
    else if (value < 0.0d0 .or. ieee_is_nan(value)) then
        rst = ieee_value(0.0d0, ieee_quiet_nan)
    else
        rst = log(value)
    end if
end function

pure function mvd_log_pdf(this, x) result(rst)
    !! Compatibility log-density fallback for custom multivariate distributions.
    !! Override this binding to compute log densities without PDF underflow.
    class(multivariate_distribution), intent(in) :: this
        !! Multivariate distribution whose PDF is evaluated.
    real(real64), intent(in), dimension(:) :: x
        !! Multivariate evaluation point.
    real(real64) :: rst
        !! Log PDF; zero density gives -infinity, negative/NaN density gives NaN.

    rst = probability_log(this%pdf(x))
end function

pure elemental function dist_log_pdf(this, x) result(rst)
    !! Computes the natural log density or probability mass.
    !! Built-in distributions use analytic log-space expressions. Other
    !! distributions fall back to log(pdf(x)); override this binding to avoid
    !! underflow in custom PDFs. Subclasses changing a built-in law must also
    !! override its log-density and tail bindings.
    class(distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point; discrete masses require integer-valued x.
    real(real64) :: rst
        !! Log density; outside support gives -infinity, invalid inputs give NaN.
        !! An integrable density singularity at a support endpoint gives +infinity.
    real(real64) :: first, second, count, argument

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (ieee_is_nan(x)) return
    select type (this)
    class is (t_distribution)
        if (.not.ieee_is_finite(this%dof) .or. this%dof <= 0.0d0) return
        argument = 0.0d0
        if (x /= 0.0d0) argument = softplus(2.0d0 * log(abs(x)) - log(this%dof))
        rst = log_gamma(0.5d0 * this%dof + 0.5d0) - log_gamma(0.5d0 * this%dof) &
            - 0.5d0 * (log(this%dof) + log(pi)) - (0.5d0 * this%dof + 0.5d0) * argument
    class is (f_distribution)
        if (.not.ieee_is_finite(this%d1) .or. .not.ieee_is_finite(this%d2)) return
        if (this%d1 <= 0.0d0 .or. this%d2 <= 0.0d0) return
        rst = ieee_value(0.0d0, ieee_negative_inf)
        if (x < 0.0d0 .or. .not.ieee_is_finite(x)) return
        first = 0.5d0 * this%d1
        second = 0.5d0 * this%d2
        if (x == 0.0d0) then
            if (first < 1.0d0) rst = ieee_value(0.0d0, ieee_positive_inf)
            if (first == 1.0d0) rst = 0.0d0
            return
        end if
        argument = log(this%d1) - log(this%d2)
        rst = first * argument + (first - 1.0d0) * log(x) &
            - (first + second) * softplus(argument + log(x)) &
            - log_gamma(first) - log_gamma(second) + log_gamma(first + second)
    class is (chi_squared_distribution)
        if (this%dof <= 0) return
        rst = ieee_value(0.0d0, ieee_negative_inf)
        if (x < 0.0d0 .or. .not.ieee_is_finite(x)) return
        first = 0.5d0 * real(this%dof, real64)
        if (x == 0.0d0) then
            if (first < 1.0d0) rst = ieee_value(0.0d0, ieee_positive_inf)
            if (first == 1.0d0) rst = -log(2.0d0)
            return
        end if
        rst = (first - 1.0d0) * log(x) - 0.5d0 * x - first * log(2.0d0) - log_gamma(first)
    class is (binomial_distribution)
        if (this%n < 0 .or. .not.ieee_is_finite(this%p)) return
        if (this%p < 0.0d0 .or. this%p > 1.0d0) return
        rst = ieee_value(0.0d0, ieee_negative_inf)
        count = real(this%n, real64)
        if (x < 0.0d0 .or. x > count) return
        if (x /= aint(x)) return
        if (this%p == 0.0d0 .or. this%p == 1.0d0) then
            if (x == count * this%p) rst = 0.0d0
            return
        end if
        rst = log_gamma(count + 1.0d0) - log_gamma(x + 1.0d0) - log_gamma(count - x + 1.0d0) &
            + x * log(this%p) + (count - x) * log_one_plus(-this%p)
    class is (poisson_distribution)
        if (.not.ieee_is_finite(this%occrence_rate) .or. this%occrence_rate < 0.0d0) return
        rst = ieee_value(0.0d0, ieee_negative_inf)
        if (x < 0.0d0 .or. .not.ieee_is_finite(x)) return
        if (x /= aint(x)) return
        if (this%occrence_rate == 0.0d0) then
            if (x == 0.0d0) rst = 0.0d0
        else
            rst = x * log(this%occrence_rate) - this%occrence_rate - log_gamma(x + 1.0d0)
        end if
    class default
        rst = probability_log(this%pdf(x))
    end select
end function

pure elemental function dist_survival(this, x) result(rst)
    !! Computes the upper-tail probability P(X > x).
    !! Built-ins use beta/gamma upper tails rather than subtracting the CDF.
    !! Custom distributions fall back to 1 - cdf(x), which can lose tail accuracy.
    class(distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point; discrete tails include masses strictly above x.
    real(real64) :: rst
        !! Probability on [0, 1]; invalid parameters or NaN x give NaN.
    real(real64) :: argument, tail, count

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (ieee_is_nan(x)) return
    select type (this)
    class is (t_distribution)
        if (.not.ieee_is_finite(this%dof) .or. this%dof <= 0.0d0) return
        argument = 1.0d0
        if (x /= 0.0d0) argument = exp(-softplus(2.0d0 * log(abs(x)) - log(this%dof)))
        tail = 0.5d0 * regularized_beta(0.5d0 * this%dof, 0.5d0, argument)
        rst = tail
        if (x < 0.0d0) rst = 1.0d0 - tail
    class is (f_distribution)
        if (.not.ieee_is_finite(this%d1) .or. .not.ieee_is_finite(this%d2)) return
        if (this%d1 <= 0.0d0 .or. this%d2 <= 0.0d0) return
        rst = 1.0d0
        if (x <= 0.0d0) return
        argument = exp(-softplus(log(this%d1) + log(x) - log(this%d2)))
        rst = regularized_beta(0.5d0 * this%d2, 0.5d0 * this%d1, argument)
    class is (chi_squared_distribution)
        if (this%dof <= 0) return
        rst = 1.0d0
        if (x <= 0.0d0) return
        rst = regularized_gamma_upper(0.5d0 * real(this%dof, real64), 0.5d0 * x)
    class is (binomial_distribution)
        if (this%n < 0 .or. .not.ieee_is_finite(this%p)) return
        if (this%p < 0.0d0 .or. this%p > 1.0d0) return
        count = real(this%n, real64)
        if (x < 0.0d0) then
            rst = 1.0d0
        else if (x >= count) then
            rst = 0.0d0
        else
            rst = regularized_beta(aint(x) + 1.0d0, count - aint(x), this%p)
        end if
    class is (poisson_distribution)
        if (.not.ieee_is_finite(this%occrence_rate) .or. this%occrence_rate < 0.0d0) return
        if (x < 0.0d0) then
            rst = 1.0d0
        else if (.not.ieee_is_finite(x)) then
            rst = 0.0d0
        else
            rst = regularized_gamma_lower(aint(x) + 1.0d0, this%occrence_rate)
        end if
    class default
        rst = 1.0d0 - this%cdf(x)
    end select
end function

pure elemental function dist_log_cdf(this, x) result(rst)
    !! Computes log(P(X <= x)), using log beta/gamma ratios for built-in laws.
    !! Custom distributions fall back to log(cdf(x)); override for extreme tails.
    class(distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point.
    real(real64) :: rst
        !! Log probability; zero probability gives -infinity, invalid inputs NaN.
    real(real64) :: argument

    rst = probability_log(this%cdf(x))
    if (ieee_is_nan(rst)) return
    select type (this)
    class is (t_distribution)
        rst = this%log_survival(-x)
    class is (f_distribution)
        if (x > 0.0d0) then
            argument = exp(-softplus(log(this%d2) - log(this%d1) - log(x)))
            rst = log_regularized_beta(0.5d0 * this%d1, 0.5d0 * this%d2, argument)
        end if
    class is (chi_squared_distribution)
        if (x > 0.0d0) rst = log_regularized_gamma_lower(0.5d0 * real(this%dof, real64), 0.5d0 * x)
    class is (binomial_distribution)
        if (x >= 0.0d0 .and. x < real(this%n, real64)) &
            rst = log_regularized_beta(real(this%n, real64) - aint(x), aint(x) + 1.0d0, 1.0d0 - this%p)
    class is (poisson_distribution)
        if (x >= 0.0d0 .and. ieee_is_finite(x)) &
            rst = log_regularized_gamma_upper(aint(x) + 1.0d0, this%occrence_rate)
    end select
end function

pure elemental function dist_log_survival(this, x) result(rst)
    !! Computes log(P(X > x)), evaluating built-in beta/gamma tails in log space.
    !! Custom distributions fall back to log(survival(x)).
    class(distribution), intent(in) :: this
        !! Distribution with valid parameters.
    real(real64), intent(in) :: x
        !! Evaluation point.
    real(real64) :: rst
        !! Log probability; zero probability gives -infinity, invalid inputs NaN.
    real(real64) :: argument, tail

    rst = probability_log(this%survival(x))
    if (ieee_is_nan(rst)) return
    select type (this)
    class is (t_distribution)
        argument = 1.0d0
        if (x /= 0.0d0) argument = exp(-softplus(2.0d0 * log(abs(x)) - log(this%dof)))
        tail = log_regularized_beta(0.5d0 * this%dof, 0.5d0, argument) - log(2.0d0)
        rst = tail
        if (x < 0.0d0) rst = stable_log1p(-exp(tail))
    class is (f_distribution)
        if (x > 0.0d0) then
            argument = exp(-softplus(log(this%d1) + log(x) - log(this%d2)))
            rst = log_regularized_beta(0.5d0 * this%d2, 0.5d0 * this%d1, argument)
        end if
    class is (chi_squared_distribution)
        if (x > 0.0d0) rst = log_regularized_gamma_upper(0.5d0 * real(this%dof, real64), 0.5d0 * x)
    class is (binomial_distribution)
        if (x >= 0.0d0 .and. x < real(this%n, real64)) &
            rst = log_regularized_beta(aint(x) + 1.0d0, real(this%n, real64) - aint(x), this%p)
    class is (poisson_distribution)
        if (x >= 0.0d0 .and. ieee_is_finite(x)) &
            rst = log_regularized_gamma_lower(aint(x) + 1.0d0, this%occrence_rate)
    end select
end function

pure elemental function nd_log_pdf(this, x) result(rst)
    !! Computes a normal log density directly, without first evaluating its PDF.
    class(normal_distribution), intent(in) :: this
        !! Normal law with finite mean and finite, positive standard deviation.
    real(real64), intent(in) :: x
        !! Evaluation point; infinite x gives -infinity.
    real(real64) :: rst
        !! Log density; invalid parameters or NaN x give NaN.
    real(real64) :: standardized

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (ieee_is_nan(x)) return
    if (.not.ieee_is_finite(this%mean_value)) return
    if (.not.ieee_is_finite(this%standard_deviation)) return
    if (this%standard_deviation <= 0.0d0) return
    standardized = normal_standardized(x, this%mean_value, this%standard_deviation) / sqrt(2.0d0)
    if (abs(standardized) > sqrt(huge(1.0d0))) then
        rst = ieee_value(0.0d0, ieee_negative_inf)
    else
        rst = -standardized**2 - log(this%standard_deviation) - 0.5d0 * log(2.0d0 * pi)
    end if
end function

pure elemental function normal_standardized(x, mu, sigma) result(rst)
    !! Standardizes a normal evaluation point without an overflowing subtraction.
    real(real64), intent(in) :: x
        !! Evaluation point; NaN inputs propagate.
    real(real64), intent(in) :: mu
        !! Finite normal mean.
    real(real64), intent(in) :: sigma
        !! Finite, strictly positive normal standard deviation.
    real(real64) :: rst
        !! Standardized point (x - mu) / sigma.

    if ((x >= 0.0d0 .and. mu >= 0.0d0) .or. (x <= 0.0d0 .and. mu <= 0.0d0)) then
        rst = (x - mu) / sigma
    else
        rst = x / sigma - mu / sigma
    end if
end function

pure elemental function nd_survival(this, x) result(rst)
    !! Computes P(X > x) for a normal law using the complementary error function.
    class(normal_distribution), intent(in) :: this
        !! Normal law with finite mean and finite, positive standard deviation.
    real(real64), intent(in) :: x
        !! Evaluation point.
    real(real64) :: rst
        !! Upper-tail probability; invalid parameters or NaN x give NaN.

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (ieee_is_nan(x)) return
    if (.not.ieee_is_finite(this%mean_value)) return
    if (.not.ieee_is_finite(this%standard_deviation)) return
    if (this%standard_deviation <= 0.0d0) return
    rst = 0.5d0 * erfc(normal_standardized(x, this%mean_value, this%standard_deviation) / sqrt(2.0d0))
end function

pure elemental function normal_log_tail(z) result(rst)
    !! Computes the standard-normal log upper tail.
    !! Uses ERFC centrally and an asymptotic Mills-ratio expansion for z >= 26,
    !! retaining finite log probabilities after the ordinary tail underflows.
    real(real64), intent(in) :: z
        !! Standardized evaluation point.
    real(real64) :: rst
        !! Log upper-tail probability; +infinity gives -infinity, NaN propagates.
    real(real64) :: term, total, previous, inverse_square
    integer(int32) :: iteration

    if (ieee_is_nan(z)) then
        rst = z
    else if (z < 0.0d0) then
        rst = stable_log1p(-0.5d0 * erfc(-z / sqrt(2.0d0)))
    else if (z < 26.0d0) then
        rst = probability_log(0.5d0 * erfc(z / sqrt(2.0d0)))
    else if (.not.ieee_is_finite(z) .or. z / sqrt(2.0d0) > sqrt(huge(1.0d0))) then
        rst = ieee_value(0.0d0, ieee_negative_inf)
    else
        inverse_square = (1.0d0 / z)**2
        term = 1.0d0
        total = term
        do iteration = 1, 100
            previous = term
            term = -term * real(2 * iteration - 1, real64) * inverse_square
            if (abs(term) >= abs(previous)) exit
            total = total + term
            if (abs(term) <= epsilon(total) * abs(total)) exit
        end do
        rst = -(z / sqrt(2.0d0))**2 - log(z) - 0.5d0 * log(2.0d0 * pi) + log(total)
    end if
end function

pure elemental function nd_log_cdf(this, x) result(rst)
    !! Computes the normal log lower tail, including extreme negative tails.
    class(normal_distribution), intent(in) :: this
        !! Normal law with finite mean and finite, positive standard deviation.
    real(real64), intent(in) :: x
        !! Evaluation point.
    real(real64) :: rst
        !! Log CDF; invalid parameters or NaN x give NaN.

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (ieee_is_nan(x)) return
    if (.not.ieee_is_finite(this%mean_value)) return
    if (.not.ieee_is_finite(this%standard_deviation)) return
    if (this%standard_deviation <= 0.0d0) return
    rst = normal_log_tail(-normal_standardized(x, this%mean_value, this%standard_deviation))
end function

pure elemental function nd_log_survival(this, x) result(rst)
    !! Computes the normal log upper tail, including extreme positive tails.
    class(normal_distribution), intent(in) :: this
        !! Normal law with finite mean and finite, positive standard deviation.
    real(real64), intent(in) :: x
        !! Evaluation point.
    real(real64) :: rst
        !! Log survival probability; invalid parameters or NaN x give NaN.

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (ieee_is_nan(x)) return
    if (.not.ieee_is_finite(this%mean_value)) return
    if (.not.ieee_is_finite(this%standard_deviation)) return
    if (this%standard_deviation <= 0.0d0) return
    rst = normal_log_tail(normal_standardized(x, this%mean_value, this%standard_deviation))
end function

pure elemental function lnd_log_pdf(this, x) result(rst)
    !! Computes a log-normal log density from the underlying normal log density.
    class(log_normal_distribution), intent(in) :: this
        !! Law with finite log-mean and finite, positive log-standard deviation.
    real(real64), intent(in) :: x
        !! Evaluation point in observation space, not log space.
    real(real64) :: rst
        !! Log density; x <= 0 gives -infinity, invalid parameters or NaN x NaN.
    type(normal_distribution) :: normal

    normal%mean_value = this%mean_value
    normal%standard_deviation = this%standard_deviation
    rst = normal%log_pdf(0.0d0)
    if (ieee_is_nan(rst) .or. ieee_is_nan(x)) then
        rst = ieee_value(0.0d0, ieee_quiet_nan)
        return
    end if
    if (x <= 0.0d0 .or. .not.ieee_is_finite(x)) then
        rst = ieee_value(0.0d0, ieee_negative_inf)
    else
        rst = normal%log_pdf(log(x)) - log(x)
    end if
end function

pure elemental function lnd_survival(this, x) result(rst)
    !! Computes the log-normal upper-tail probability P(X > x).
    class(log_normal_distribution), intent(in) :: this
        !! Law with finite log-mean and finite, positive log-standard deviation.
    real(real64), intent(in) :: x
        !! Evaluation point in observation space.
    real(real64) :: rst
        !! Survival probability; x <= 0 gives one, invalid inputs give NaN.

    rst = exp(this%log_survival(x))
end function

pure elemental function lnd_log_cdf(this, x) result(rst)
    !! Computes the log-normal log lower tail using the normal log-CDF.
    class(log_normal_distribution), intent(in) :: this
        !! Law with finite log-mean and finite, positive log-standard deviation.
    real(real64), intent(in) :: x
        !! Evaluation point in observation space.
    real(real64) :: rst
        !! Log CDF; x <= 0 gives -infinity, invalid inputs give NaN.
    type(normal_distribution) :: normal

    normal%mean_value = this%mean_value
    normal%standard_deviation = this%standard_deviation
    rst = normal%log_cdf(0.0d0)
    if (ieee_is_nan(rst) .or. ieee_is_nan(x)) then
        rst = ieee_value(0.0d0, ieee_quiet_nan)
    else if (x <= 0.0d0) then
        rst = ieee_value(0.0d0, ieee_negative_inf)
    else
        rst = normal%log_cdf(log(x))
    end if
end function

pure elemental function lnd_log_survival(this, x) result(rst)
    !! Computes the log-normal log upper tail using the normal log-survival.
    class(log_normal_distribution), intent(in) :: this
        !! Law with finite log-mean and finite, positive log-standard deviation.
    real(real64), intent(in) :: x
        !! Evaluation point in observation space.
    real(real64) :: rst
        !! Log survival probability; x <= 0 gives zero, invalid inputs give NaN.
    type(normal_distribution) :: normal

    normal%mean_value = this%mean_value
    normal%standard_deviation = this%standard_deviation
    rst = normal%log_survival(0.0d0)
    if (ieee_is_nan(rst) .or. ieee_is_nan(x)) then
        rst = ieee_value(0.0d0, ieee_quiet_nan)
    else if (x <= 0.0d0) then
        rst = 0.0d0
    else
        rst = normal%log_survival(log(x))
    end if
end function

! ------------------------------------------------------------------------------
pure elemental function dist_std_var(this, x) result(rst)
    !! Computes the standardized variable for the distribution.
    class(distribution), intent(in) :: this
        !! The distribution object.
    real(real64), intent(in) :: x
        !! The value of interest.
    real(real64) :: rst
        !! The result.

    ! Local Variables
    integer(int32), parameter :: maxiter = 100
    real(real64), parameter :: tol = 1.0d-6
    integer(int32) :: i
    real(real64) :: f, df, h, twoh, dy

    ! Process
    !
    ! We use a simplified Newton's method to solve for the independent variable
    ! of the CDF function
    h = 1.0d-6
    twoh = 2.0d0 * h
    rst = 0.5d0 ! just an initial guess
    do i = 1, maxiter
        ! Compute the CDF and its derivative at y
        f = this%cdf(rst) - x
        df = (this%cdf(rst + h) - this%cdf(rst - h)) / twoh
        dy = f / df
        rst = rst - dy
        if (abs(dy) < tol) exit
    end do
end function

! ------------------------------------------------------------------------------
pure function dist_defined_range(this) result(rst)
    !! Gets the defined range for the distribution.
    class(distribution), intent(in) :: this
        !! The distribution object.
    real(real64), dimension(2) :: rst
        !! The defined range of the probability distributions [min, max].  In
        !! the event that either min or max are infinite, a value of huge(0.0d0)
        !! is returned as opposed to infinite to avoid possible issues with
        !! using these values directly.

    rst = [-huge(0.0d0), huge(0.0d0)]
end function

! ******************************************************************************
! STUDENT'S T-DISTRIBUTION
! ------------------------------------------------------------------------------
! REF: https://en.wikipedia.org/wiki/Student%27s_t-distribution
pure elemental function td_pdf(this, x) result(rst)
    !! Computes the probability density function.
    !!
    !! The PDF for Student's T-Distribution is given as 
    !! $$ f(t) = \frac{ \Gamma \left( \frac{\nu + 1}{2} \right) }
    !! { \sqrt{\nu \pi} \Gamma \left( \frac{\nu}{2} \right) } 
    !! \left( 1 + \frac{t^2}{\nu} \right)^{-(\nu + 1) / 2} $$.
    class(t_distribution), intent(in) :: this
        !! The t_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    ! Process
    rst = exp(this%log_pdf(x))
end function

! ------------------------------------------------------------------------------
pure elemental function td_cdf(this, x) result(rst)
    !! Computes the cumulative distribution function.
    !!
    !! The CDF for Student's T-Distribution is given as
    !! $$ F(t) = \int_{-\infty}^{t} f(u) \,du = 1 - \frac{1}{2} I_{x(t)}
    !! \left( \frac{\nu}{2}, \frac{1}{2} \right) $$
    !! where $$ x(t) = \frac{\nu}{\nu + t^2} $$.
    class(t_distribution), intent(in) :: this
        !! The t_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    ! Process
    real(real64) :: t
    rst = this%survival(-x)
end function

! ------------------------------------------------------------------------------
pure function td_mean(this) result(rst)
    !! Computes the mean of the distribution.
    class(t_distribution), intent(in) :: this
        !! The t_distribution object.
    real(real64) :: rst
        !! The mean.

    ! Process
    if (this%dof < 1.0d0) then
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    else
        rst = 0.0d0
    end if
end function

! ------------------------------------------------------------------------------
pure function td_median(this) result(rst)
    !! Computes the median of the distribution.
    class(t_distribution), intent(in) :: this
        !! The t_distribution object.
    real(real64) :: rst

    ! Process
    rst = 0.0d0
end function

! ------------------------------------------------------------------------------
pure function td_mode(this) result(rst)
    !! Computes the mode of the distribution.
    class(t_distribution), intent(in) :: this
        !! The t_distribution object.
    real(real64) :: rst
        !! The mode.

    ! Process
    rst = 0.0d0
end function

! ------------------------------------------------------------------------------
pure function td_variance(this) result(rst)
    !! Computes the variance of the distribution.
    class(t_distribution), intent(in) :: this
        !! The t_distribution object.
    real(real64) :: rst
        !! The variance.

    ! Process
    if (this%dof <= 1.0d0) then
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    else if (this%dof > 1.0d0 .and. this%dof <= 2.0d0) then
        rst = ieee_value(rst, IEEE_POSITIVE_INF)
    else
        rst = this%dof / (this%dof - 2.0d0)
    end if
end function

! ------------------------------------------------------------------------------
subroutine td_recenter(this, x)
    !! Recenters the distribution about the supplied value.  This routine has
    !! no effect for this distribution as it is always centered about 0.
    class(t_distribution), intent(inout) :: this
        !! The t_distribution object.
    real(real64), intent(in) :: x
        !! The value about which to recenter.
end subroutine

! ******************************************************************************
! NORMAL DISTRIBUTION
! ------------------------------------------------------------------------------
pure elemental function nd_pdf(this, x) result(rst)
    !! Computes the probability density function.
    !!
    !! The PDF for a normal distribution is given as 
    !! $$ f(x) = \frac{1}{\sigma \sqrt{2 \pi}} \exp \left(-\frac{1}{2} 
    !! \left( \frac{x - \mu}{\sigma} \right)^2 \right) $$.
    class(normal_distribution), intent(in) :: this
        !! The normal_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    rst = exp(this%log_pdf(x))
end function

! ------------------------------------------------------------------------------
pure elemental function nd_cdf(this, x) result(rst)
    !! Computes the cumulative distribution function.
    !!
    !! The CDF for a normal distribution is given as 
    !! $$ F(x) = \frac{1}{2} \left( 1 + erf \left( \frac{x - \mu}
    !! {\sigma \sqrt{2}} \right) \right) $$.
    class(normal_distribution), intent(in) :: this
        !! The normal_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    rst = exp(this%log_cdf(x))
end function

! ------------------------------------------------------------------------------
pure function nd_mean(this) result(rst)
    !! Computes the mean of the distribution.
    class(normal_distribution), intent(in) :: this
        !! The normal_distribution object.
    real(real64) :: rst
        !! The mean
    rst = this%mean_value
end function

! ------------------------------------------------------------------------------
pure function nd_median(this) result(rst)
    !! Computes the median of the distribution.
    class(normal_distribution), intent(in) :: this
        !! The normal_distribution object.
    real(real64) :: rst
        !! The median.
    rst = this%mean_value
end function

! ------------------------------------------------------------------------------
pure function nd_mode(this) result(rst)
    !! Computes the mode of the distribution.
    class(normal_distribution), intent(in) :: this
        !! The normal_distribution object.
    real(real64) :: rst
        !! The mode.
    rst = this%mean_value
end function

! ------------------------------------------------------------------------------
pure function nd_variance(this) result(rst)
    !! Computes the variance of the distribution.
    class(normal_distribution), intent(in) :: this
        !! The normal_distribution object.
    real(real64) :: rst
        !! The variance.
    rst = this%standard_deviation**2
end function

! ------------------------------------------------------------------------------
subroutine nd_standardize(this)
    !! Standardizes the normal distribution to a mean of 0 and a 
    !! standard deviation of 1.
    class(normal_distribution), intent(inout) :: this
        !! The normal_distribution object.
    this%mean_value = 0.0d0
    this%standard_deviation = 1.0d0
end subroutine

! ------------------------------------------------------------------------------
subroutine nd_recenter(this, x)
    !! Recenters the distribution about the supplied value.
    class(normal_distribution), intent(inout) :: this
        !! The normal_distribution object.
    real(real64), intent(in) :: x
        !! The value about which to recenter.

    ! Process
    this%mean_value = x
end subroutine

! ******************************************************************************
! F DISTRIBUTION
! ------------------------------------------------------------------------------
pure elemental function fd_pdf(this, x) result(rst)
    !! Computes the probability density function.
    !!
    !! The PDF for a F distribution is given as 
    !! $$ f(x) = 
    !! \sqrt{ \frac{ (d_1 x)^{d_1} d_{2}^{d_2} }{ (d_1 x + d_2)^{d_1 + d_2} } } 
    !! \frac{1}{x \beta \left( \frac{d_1}{2}, \frac{d_2}{2} \right) } $$.
    class(f_distribution), intent(in) :: this
        !! The f_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    ! Process
    real(real64) :: d1, d2
    d1 = this%d1
    d2 = this%d2
    rst = exp(this%log_pdf(x))
end function

! ------------------------------------------------------------------------------
pure elemental function fd_cdf(this, x) result(rst)
    !! Computes the cumulative distribution function.
    !!
    !! The CDF for a F distribution is given as 
    !! $$ F(x) = I_{d_1 x/(d_1 x + d_2)} \left( \frac{d_1}{2}, 
    !! \frac{d_2}{2} \right) $$.
    class(f_distribution), intent(in) :: this
        !! The f_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    ! Process
    real(real64) :: d1, d2
    d1 = this%d1
    d2 = this%d2
    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (.not.ieee_is_finite(d1) .or. .not.ieee_is_finite(d2)) return
    if (d1 <= 0.0d0 .or. d2 <= 0.0d0 .or. ieee_is_nan(x)) return
    if (x <= 0.0d0) then
        rst = 0.0d0
    else
        rst = regularized_beta(0.5d0 * d1, 0.5d0 * d2, &
            exp(-softplus(log(d2) - log(d1) - log(x))))
    end if
end function

! ------------------------------------------------------------------------------
pure function fd_mean(this) result(rst)
    !! Computes the mean of the distribution.
    class(f_distribution), intent(in) :: this
        !! The f_distribution object.
    real(real64) :: rst
        !! The mean.

    ! Process
    if (this%d2 > 2.0d0) then
        rst = this%d2 / (this%d2 - 2.0d0)
    else
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    end if
end function

! ------------------------------------------------------------------------------
pure function fd_median(this) result(rst)
    !! Computes the median of the distribution.
    class(f_distribution), intent(in) :: this
        !! The f_distribution object.
    real(real64) :: rst
        !! The median.
    rst = ieee_value(rst, IEEE_QUIET_NAN)
end function

! ------------------------------------------------------------------------------
pure function fd_mode(this) result(rst)
    !! Computes the mode of the distribution.
    class(f_distribution), intent(in) :: this
        !! The f_distribution object.
    real(real64) :: rst
        !! The mode.

    ! Process
    if (this%d1 > 2.0d0) then
        rst = ((this%d1 - 2.0d0) / this%d1) * (this%d2 / (this%d2 + 2.0d0))
    else
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    end if
end function

! ------------------------------------------------------------------------------
pure function fd_variance(this) result(rst)
    !! Computes the variance of the distribution.
    class(f_distribution), intent(in) :: this
        !! The f_distribution object.
    real(real64) :: rst
        !! The variance.

    ! Process
    real(real64) :: d1, d2
    d1 = this%d1
    d2 = this%d2
    if (d2 > 4.0d0) then
        rst = (2.0d0 * d2**2 * (d1 + d2 - 2.0d0)) / &
            (d1 * (d2 - 2.0d0)**2 * (d2 - 4.0d0))
    else
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    end if
end function

! ------------------------------------------------------------------------------
pure function fd_range(this) result(rst)
    !! Gets the defined range for the distribution.
    class(f_distribution), intent(in) :: this
        !! The f_distribution object.
    real(real64), dimension(2) :: rst
        !! The defined range of the probability distributions [0, infinity).  As
        !! using a value of infinity may cause issue, this routine returns
        !! huge(0.0d0) instead.

    rst = [0.0d0, huge(0.0d0)]
end function

! ------------------------------------------------------------------------------
subroutine fd_recenter(this, x)
    !! Recenters the distribution about the supplied value.
    class(f_distribution), intent(inout) :: this
        !! The f_distribution object.
    real(real64), intent(in) :: x
        !! The value about which to recenter.

    ! Process
    this%d2 = 2.0d0 * x / (x - 1.0d0)
end subroutine

! ******************************************************************************
! CHI-SQUARED DISTRIBUTION
! ------------------------------------------------------------------------------
pure elemental function cs_pdf(this, x) result(rst)
    !! Computes the probability density function.
    !!
    !! The PDF for a Chi-squared distribution is given as 
    !! $$ f(x) = \frac{x^{k/2 - 1} \exp{-x / 2}} {2^{k / 2} 
    !! \Gamma \left( \frac{k}{2} \right)} $$.
    class(chi_squared_distribution), intent(in) :: this
        !! The chi_squared_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    ! Local Variables
    real(real64) :: arg

    ! Process
    arg = 0.5d0 * this%dof
    rst = exp(this%log_pdf(x))
end function

! ------------------------------------------------------------------------------
pure elemental function cs_cdf(this, x) result(rst)
    !! Computes the cumulative distribution function.
    !!
    !! The CDF for a Chi-squared distribution is given as 
    !! $$ F(x) = \frac{ \gamma \left( \frac{k}{2}, \frac{x}{2} \right) }
    !! { \Gamma \left( \frac{k}{2} \right)} $$.
    class(chi_squared_distribution), intent(in) :: this
        !! The chi_squared_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    ! Local Variables
    real(real64) :: arg

    ! Process
    arg = 0.5d0 * this%dof
    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (this%dof <= 0 .or. ieee_is_nan(x)) return
    if (x <= 0.0d0) then
        rst = 0.0d0
    else
        rst = regularized_gamma_lower(arg, 0.5d0 * x)
    end if
end function

! ------------------------------------------------------------------------------
pure function cs_mean(this) result(rst)
    !! Computes the mean of the distribution.
    class(chi_squared_distribution), intent(in) :: this
        !! The chi_squared_distribution object.
    real(real64) :: rst
        !! The mean.

    ! Process
    rst = real(this%dof, real64)
end function

! ------------------------------------------------------------------------------
pure function cs_median(this) result(rst)
    !! Computes the median of the distribution.
    class(chi_squared_distribution), intent(in) :: this
        !! The chi_squared_distribution object.
    real(real64) :: rst
        !! The median.

    ! Process
    rst = this%dof * (1.0d0 - 2.0d0 / (9.0d0 * this%dof))**3
end function

! ------------------------------------------------------------------------------
pure function cs_mode(this) result(rst)
    !! Computes the mode of the distribution.
    class(chi_squared_distribution), intent(in) :: this
        !! The chi_squared_distribution object.
    real(real64) :: rst
        !! The mode.

    ! Process
    rst = max(this%dof - 2.0d0, 0.0d0)
end function

! ------------------------------------------------------------------------------
pure function cs_variance(this) result(rst)
    !! Computes the variance of the distribution.
    class(chi_squared_distribution), intent(in) :: this
        !! The chi_squared_distribution object.
    real(real64) :: rst
        !! The variance.

    ! Process
    rst = 2.0d0 * this%dof
end function

! ------------------------------------------------------------------------------
pure function cs_range(this) result(rst)
    !! Gets the defined range for the distribution.
    class(chi_squared_distribution), intent(in) :: this
        !! The chi_squared_distribution object.
    real(real64), dimension(2) :: rst
        !! The defined range of the probability distributions [0, infinity).  As
        !! using a value of infinity may cause issue, this routine returns
        !! huge(0.0d0) instead.

    rst = [0.0d0, huge(0.0d0)]
end function

! ------------------------------------------------------------------------------
subroutine cs_recenter(this, x)
    !! Recenters the distribution about the supplied value.
    class(chi_squared_distribution), intent(inout) :: this
        !! The chi_squared_distribution object.
    real(real64), intent(in) :: x
        !! The value about which to recenter.

    ! Process
    this%dof = floor(x)
end subroutine

! ******************************************************************************
! BINOMIAL DISTRIBUTION
! ------------------------------------------------------------------------------
pure elemental function bd_pdf(this, x) result(rst)
    !! Computes the probability mass function.
    !!
    !! The PMF for a binomial distribution is given as 
    !! $$ f(k,n,p) = \frac{n!}{k! \left( n - k! \right)} p^k 
    !! \left( 1 - p \right)^{n-k} $$.
    class(binomial_distribution), intent(in) :: this
        !! The binomial_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.  This parameter
        !! is the number k successes in the n independent trials.  As
        !! such, this parameter must exist on the set [0, n].
    real(real64) :: rst
        !! The value of the function.

    ! Local Variables
    real(real64) :: dn

    ! Process
    dn = real(this%n, real64)
    rst = exp(this%log_pdf(x))
end function

! ------------------------------------------------------------------------------
pure elemental function bd_cdf(this, x) result(rst)
    !! Computes the cumulative distribution funtion.
    !!
    !! The CDF for a binomial distribution is given as 
    !! $$ F(k,n,p) = I_{1-p} \left( n - k, 1 + k \right) $$, which is simply
    !! the regularized incomplete beta function.
    class(binomial_distribution), intent(in) :: this
        !! The binomial_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.  This parameter
        !! is the number k successes in the n independent trials.  As
        !! such, this parameter must exist on the set [0, n].
    real(real64) :: rst
        !! The value of the function.

    ! Local Variables
    real(real64) :: dn

    ! Process
    dn = real(this%n, real64)
    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (this%n < 0 .or. .not.ieee_is_finite(this%p) .or. ieee_is_nan(x)) return
    if (this%p < 0.0d0 .or. this%p > 1.0d0) return
    if (x < 0.0d0) then
        rst = 0.0d0
    else if (x >= dn) then
        rst = 1.0d0
    else
        rst = regularized_beta(dn - aint(x), aint(x) + 1.0d0, 1.0d0 - this%p)
    end if
end function

! ------------------------------------------------------------------------------
pure function bd_mean(this) result(rst)
    !! Computes the mean of the distribution.
    class(binomial_distribution), intent(in) :: this
        !! The binomial_distribution object.
    real(real64) :: rst
        !! The mean.

    rst = real(this%n * this%p, real64)
end function

! ------------------------------------------------------------------------------
pure function bd_median(this) result(rst)
    !! Computes the median of the distribution.
    class(binomial_distribution), intent(in) :: this
        !! The binomial_distribution object.
    real(real64) :: rst
        !! The median.

    rst = real(this%n * this%p, real64)
end function

! ------------------------------------------------------------------------------
pure function bd_mode(this) result(rst)
    !! Computes the mode of the distribution.
    class(binomial_distribution), intent(in) :: this
        !! The binomial_distribution object.
    real(real64) :: rst
        !! The mode.

    rst = (this%n + 1.0d0) * this%p
end function

! ------------------------------------------------------------------------------
pure function bd_variance(this) result(rst)
    !! Computes the variance of the distribution.
    class(binomial_distribution), intent(in) :: this
        !! The binomial_distribution object.
    real(real64) :: rst
        !! The variance.

    rst = this%n * this%p * (1.0d0 - this%p)
end function

! ------------------------------------------------------------------------------
pure function bd_range(this) result(rst)
    !! Gets the defined range for the distribution.
    class(binomial_distribution), intent(in) :: this
        !! The binomial_distribution object.
    real(real64), dimension(2) :: rst
        !! The defined range of the probability distributions [0, infinity).  As
        !! using a value of infinity may cause issue, this routine returns
        !! huge(0.0d0) instead.

    rst = [0.0d0, huge(0.0d0)]
end function

! ------------------------------------------------------------------------------
subroutine bd_recenter(this, x)
    !! Recenters the distribution about the supplied value.
    class(binomial_distribution), intent(inout) :: this
        !! The binomial_distribution object.
    real(real64), intent(in) :: x
        !! The value about which to recenter.

    ! Process
    this%p = x / this%n
end subroutine

! ******************************************************************************
! MULTIVARIATE NORMAL DISTRIBUTION
! ------------------------------------------------------------------------------
pure subroutine mvnd_init(this, mu, sigma)
    use linalg, only : cholesky_factor
    !! Initializes the multivariate normal distribution by defining the mean
    !! values and covariance matrix.
    class(multivariate_normal_distribution), intent(inout) :: this
        !! The multivariate_normal_distribution object.
    real(real64), intent(in), dimension(:) :: mu
        !! An N-element array containing the mean values.
    real(real64), intent(in), dimension(:,:) :: sigma
        !! The N-by-N covariance matrix.  The PDF exists only if this matrix
        !! is positive-definite; therefore, the positive-definite constraint 
        !! is checked within this routine and enforced.  An error is thrown if
        !! the supplied matrix is not positive-definite.

    ! Local Variables
    integer(int32) :: n, row_index
    real(real64) :: scale
    
    ! Initialization
    n = size(mu)

    ! Input Checking
    if (size(sigma, 1) /= n .or. size(sigma, 2) /= n) error stop FS_MATRIX_SIZE_ERROR
    if (n < 1) error stop FS_INVALID_INPUT_ERROR
    if (.not.all(ieee_is_finite(mu)) .or. .not.all(ieee_is_finite(sigma))) &
        error stop FS_INVALID_INPUT_ERROR
    scale = maxval(abs(sigma))
    if (scale == 0.0d0) error stop FS_INVALID_INPUT_ERROR
    if (any(abs(sigma / scale - transpose(sigma) / scale) > &
        10.0d0 * real(n, real64) * epsilon(scale))) error stop FS_INVALID_INPUT_ERROR

    ! Store the matrices
    this%m_means = mu
    this%m_cov = sigma
    ! Compute the Cholesky factorization of the covariance matrix
    this%m_cholesky = cholesky_factor(sigma, upper = .false.)
    if (.not.all(ieee_is_finite(this%m_cholesky))) error stop FS_INVALID_INPUT_ERROR
    this%m_logCovDet = 0.0d0
    do row_index = 1, n
        this%m_logCovDet = this%m_logCovDet + 2.0d0 * log(this%m_cholesky(row_index, row_index))
    end do
end subroutine

! ------------------------------------------------------------------------------
pure function mvnd_pdf(this, x) result(rst)
    !! Evaluates the PDF for the multivariate normal distribution.
    class(multivariate_normal_distribution), intent(in) :: this
        !! The multivariate_normal_distribution object.
    real(real64), intent(in), dimension(:) :: x
        !! The values at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    rst = exp(this%log_pdf(x))
end function

pure function mvnd_log_pdf(this, x) result(rst)
    !! Computes the multivariate-normal log density without an explicit inverse.
    !! A triangular solve whitens the residual; normalization uses a log determinant.
    use blas, only : dtrsv
    class(multivariate_normal_distribution), intent(in) :: this
        !! Initialized distribution with a symmetric positive-definite covariance.
    real(real64), intent(in), dimension(:) :: x
        !! Evaluation point, the same length as the stored mean vector.
    real(real64) :: rst
        !! Log density; uninitialized objects or nonfinite points give NaN.
        !! An unrepresentably large squared residual gives negative infinity.
    real(real64), allocatable, dimension(:) :: delta
    real(real64) :: residual_norm
    integer(int32) :: n

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (.not.allocated(this%m_means)) return
    n = size(this%m_means)
    if (size(x) /= n) error stop FS_ARRAY_SIZE_ERROR
    if (.not.all(ieee_is_finite(x))) return
    delta = x - this%m_means
    call dtrsv('L', 'N', 'N', n, this%m_cholesky, n, delta, 1)
    residual_norm = norm2(delta) / sqrt(2.0d0)
    if (residual_norm > sqrt(huge(1.0d0))) then
        rst = ieee_value(0.0d0, ieee_negative_inf)
    else
        rst = -residual_norm**2 - 0.5d0 * &
            (real(n, real64) * log(2.0d0 * pi) + this%m_logCovDet)
    end if
end function

! ------------------------------------------------------------------------------
pure subroutine mvnd_update_mean(this, x)
    !! Updates the mean value array.
    class(multivariate_normal_distribution), intent(inout) :: this
        !! The multivariate_normal_distribution object.
    real(real64), intent(in), dimension(:) :: x
        !! The N-element array of new mean values.

    ! Local Variables
    integer(int32) :: n, nc
    
    ! Initialization
    n = size(x)
    nc = size(this%m_means)

    ! Process
    if (.not.allocated(this%m_means)) then
        ! This is an initial set-up - just store the values and be done
        allocate(this%m_means(n), source = x)
    end if

    ! Else, ensure the array is of the correct size before updating
    if (n /= nc) error stop 2
    this%m_means = x
end subroutine

! ------------------------------------------------------------------------------
pure function mvnd_get_means(this) result(rst)
    !! Gets the mean values of the distribution.
    class(multivariate_normal_distribution), intent(in) :: this
        !! The multivariate_normal_distribution object.
    real(real64), allocatable, dimension(:) :: rst
        !! The mean values.

    ! Process
    integer(int32) :: n
    if (allocated(this%m_means)) then
        n = size(this%m_means)
        allocate(rst(n), source = this%m_means)
    else
        allocate(rst(0))
    end if
end function

! ------------------------------------------------------------------------------
pure function mvnd_get_covariance(this) result(rst)
    !! Gets the covariance matrix of the distribution.
    class(multivariate_normal_distribution), intent(in) :: this
        !! The multivariate_normal_distribution object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The covariance matrix.

    ! Process
    integer(int32) :: n
    if (allocated(this%m_cov)) then
        n = size(this%m_cov, 1)
        allocate(rst(n, n), source = this%m_cov)
    else
        allocate(rst(0, 0))
    end if
end function

! ------------------------------------------------------------------------------
pure function mvnd_get_cholesky(this) result(rst)
    !! Gets the lower triangular form of the Cholesky factorization of the
    !! covariance matrix of the distribution.
    class(multivariate_normal_distribution), intent(in) :: this
        !! The multivariate_normal_distribution object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The Cholesky factored matrix.

    ! Process
    integer(int32) :: n
    if (allocated(this%m_cholesky)) then
        n = size(this%m_cholesky, 1)
        allocate(rst(n, n), source = this%m_cholesky)
    else
        allocate(rst(0, 0))
    end if
end function

! ******************************************************************************
! LOG NORMAL DISTRIBUTION
! ------------------------------------------------------------------------------
pure elemental function lnd_pdf(this, x) result(rst)
    !! Computes the probability density function.
    !!
    !! The PDF for a log-normal distribution is given as
    !! $$ f(x) = \frac{1}{x \sigma \sqrt{2 \pi}} \exp{\left(- \frac{\left( 
    !! \ln{x} - \mu \right)^2}{2 \sigma^2} \right)} $$
    class(log_normal_distribution), intent(in) :: this
        !! The log_normal_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    rst = exp(this%log_pdf(x))
end function

! ------------------------------------------------------------------------------
pure elemental function lnd_cdf(this, x) result(rst)
    !! Computes the cumulative distribution function.
    !!
    !! The CDF for a log-normal distribution is given as
    !! $$ F(x) = \frac{1}{2} \left(1 + erf\left( \frac{\ln{x} - \mu}
    !! {\sigma \sqrt{2}} \right) \right) $$
    class(log_normal_distribution), intent(in) :: this
        !! The log_normal_distribution object.
    real(real64), intent(in) :: x
        !! The value at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    rst = exp(this%log_cdf(x))
end function

! ------------------------------------------------------------------------------
pure function lnd_mean(this) result(rst)
    !! Computes the mean of the distribution
    class(log_normal_distribution), intent(in) :: this
        !! The log_normal distribution object.
    real(real64) :: rst
        !! The mean
    rst = exp(this%mean_value + 0.5d0 * this%standard_deviation**2)
end function

! ------------------------------------------------------------------------------
pure function lnd_median(this) result(rst)
    !! Computes the median of the distribution
    class(log_normal_distribution), intent(in) :: this
        !! The log_normal distribution object.
    real(real64) :: rst
        !! The median
    rst = exp(this%mean_value)
end function

! ------------------------------------------------------------------------------
pure function lnd_mode(this) result(rst)
    !! Computes the mode of the distribution
    class(log_normal_distribution), intent(in) :: this
        !! The log_normal distribution object.
    real(real64) :: rst
        !! The mode
    rst = exp(this%mean_value - this%standard_deviation**2)
end function

! ------------------------------------------------------------------------------
pure function lnd_variance(this) result(rst)
    !! Computes the variance of the distribution
    class(log_normal_distribution), intent(in) :: this
        !! The log_normal distribution object.
    real(real64) :: rst
        !! The variance
    rst = stable_expm1(this%standard_deviation**2) * &
        exp(2.0d0 * this%mean_value + this%standard_deviation**2)
end function

! ------------------------------------------------------------------------------
pure function lnd_range(this) result(rst)
    !! Gets the defined range for the distribution.
    class(log_normal_distribution), intent(in) :: this
        !! The log_normal_distribution object.
    real(real64), dimension(2) :: rst
        !! The defined range of the probability distributions [0, infinity).  As
        !! using a value of infinity may cause issue, this routine returns
        !! huge(0.0d0) instead.

    rst = [0.0d0, huge(0.0d0)]
end function

! ------------------------------------------------------------------------------
subroutine lnd_recenter(this, x)
    !! Recenters the distribution about the supplied value.
    class(log_normal_distribution), intent(inout) :: this
        !! The log_normal_distribution object.
    real(real64), intent(in) :: x
        !! The value about which to recenter.

    this%mean_value = x
end subroutine

! ******************************************************************************
! POISSON DISTRIBUTION
! ------------------------------------------------------------------------------
pure elemental function pd_pdf(this, x) result(rst)
    !! Computes the probability mass function.
    !!
    !! The PMF for the Poisson distribution is given as \( f(k) = 
    !! \frac{\lambda^{k} e^{-\lambda}}{k!} \).
    class(poisson_distribution), intent(in) :: this
        !! The poisson_distribution object.
    real(real64), intent(in) :: x
        !! The number of occurrences (\(k\)).
    real(real64) :: rst
        !! The value of the function.
    
    ! Local Variables
    real(real64) :: lambda

    ! Process
    lambda = this%occrence_rate
    rst = exp(this%log_pdf(x))
end function

! ------------------------------------------------------------------------------
pure elemental function pd_cdf(this, x) result(rst)
    !! Computes the cumulative distribution function.
    !!
    !! The CDF for the Poisson distribution is given as 
    !! $$ F(k) = \frac{\Gamma(floor(k + 1), \lambda)}{floor(\lambda}!)} $$.
    class(poisson_distribution), intent(in) :: this
        !! The poisson_distribution object.
    real(real64), intent(in) :: x
        !! The number of occurrences (\(k\)).
    real(real64) :: rst
        !! The value of the function.

    ! Local Variables
    real(real64) :: lambda

    ! Process
    lambda = this%occrence_rate
    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (.not.ieee_is_finite(lambda) .or. ieee_is_nan(x)) return
    if (lambda < 0.0d0) return
    if (x < 0.0d0) then
        rst = 0.0d0
    else if (.not.ieee_is_finite(x)) then
        rst = 1.0d0
    else
        rst = regularized_gamma_upper(aint(x) + 1.0d0, lambda)
    end if
end function

! ------------------------------------------------------------------------------
pure function pd_mean(this) result(rst)
    !! Computes the mean of the distribution.
    class(poisson_distribution), intent(in) :: this
        !! The poisson_distribution object.
    real(real64) :: rst
        !! The mean.

    ! Process
    rst = this%occrence_rate
end function

! ------------------------------------------------------------------------------
pure function pd_median(this) result(rst)
    !! Computes the median of the distribution.
    class(poisson_distribution), intent(in) :: this
        !! The poisson_distribution object.
    real(real64) :: rst
        !! The median.

    real(real64) :: lambda
    lambda = this%occrence_rate
    rst = floor(lambda + 1.0d0 / 3.0d0 - 1.0d0 / (5.0d1 * lambda))
end function

! ------------------------------------------------------------------------------
pure function pd_mode(this) result(rst)
    !! Computes the mode of the distribution.
    class(poisson_distribution), intent(in) :: this
    !! The poisson_distribution object.
    real(real64) :: rst
        !! The mode.

    rst = max( &
        real(ceiling(this%occrence_rate) - 1.0d0, real64), &
        real(floor(this%occrence_rate), real64))
end function

! ------------------------------------------------------------------------------
pure function pd_variance(this) result(rst)
    !! Computes the variance.
    class(poisson_distribution), intent(in) :: this
        !! The poisson_distribution object.
    real(real64) :: rst
        !! The variance.

    rst = this%occrence_rate
end function

! ------------------------------------------------------------------------------
subroutine pd_recenter(this, x)
    !! Recenters the distribution about the supplied value.  This routine has
    !! no effect for this distribution as it is always centered about 0.
    class(poisson_distribution), intent(inout) :: this
        !! The poisson_distribution object.
    real(real64), intent(in) :: x
        !! The value about which to recenter.

    this%occrence_rate = x
end subroutine

! ------------------------------------------------------------------------------
pure function pd_range(this) result(rst)
    !! Gets the defined range for the distribution.
    class(poisson_distribution), intent(in) :: this
        !! The poisson_distribution object.
    real(real64), dimension(2) :: rst
        !! The defined range of the probability distributions [0, infinity).  As
        !! using a value of infinity may cause issue, this routine returns
        !! huge(0.0d0) instead.

    rst = [0.0d0, huge(0.0d0)]
end function

! ******************************************************************************
! SUPPORTING ROUTINES
! ------------------------------------------------------------------------------
pure subroutine cholesky_inverse(x, u)
    use linalg, only : solve_triangular_system
    !! Computes the inverse of a Cholesky-factored matrix.
    real(real64), intent(in), dimension(:,:) :: x
        !! The lower-triangular Cholesky factored matrix.
    real(real64), intent(inout), dimension(:,:) :: u
        !! On input, an N-by-N identity matrix.  On output, the N-by-N inverted
        !! matrix.

    ! To compute the inverse of a Cholesky factored matrix (L) consider the 
    ! following:
    !
    ! A = L * L**T
    !
    ! (L * L**T) * inv(A) = I, where I is an identity matrix
    !
    ! First, solve L * U = I, for the N-by-N matrix U
    !
    ! And then solve L' * inv(A) = U for inv(A)

    ! Solve L * U = I for U
    u = solve_triangular_system(.true., .false., .false., .true., 1.0d0, x, u)

    ! Solve L**T * inv(A) = U for inv(A)
    u = solve_triangular_system(.true., .false., .true., .true., 1.0d0, x, u)
end subroutine

! ------------------------------------------------------------------------------
pure function cholesky_determinant(x) result(rst)
    !! Computes the determinant of a Cholesky factored (lower) matrix.
    real(real64), intent(in), dimension(:,:) :: x
        !! The lower-triangular Cholesky-factored matrix.
    real(real64) :: rst
        !! The determinant.

    ! Local Variables
    integer(int32) :: i, ep, n
    real(real64) :: temp

    ! Initialization
    n = size(x, 1)
    rst = 0.0d0

    ! Compute the product of the squares of the diagonal
    temp = 1.0d0
    ep = 0
    do i = 1, n
        temp = (x(i,i))**2 * temp
        if (temp == 0.0d0) then
            rst = 0.0d0
            return
        end if

        do while (abs(temp) < 1.0d0)
            temp = 1.0d1 * temp
            ep = ep - 1
        end do

        do while (abs(temp) > 1.0d1)
            temp = 1.0d-1 * temp
            ep = ep + 1
        end do
    end do
    rst = temp * (1.0d1)**ep
end function

! ------------------------------------------------------------------------------
pure subroutine populate_identity(x)
    !! Populates the supplied matrix as an identity matrix.
    real(real64), intent(inout), dimension(:,:) :: x

    ! Local Variables
    integer(int32) :: i, m, n, mn

    ! Process
    m = size(x, 1)
    n = size(x, 2)
    mn = min(m, n)
    x = 0.0d0
    do i = 1, mn
        x(i,i) = 1.0d0
    end do
end subroutine

! ------------------------------------------------------------------------------
end module
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
    implicit none
    private
    public :: distribution
    public :: probability_log, softplus, distribution_pi
    public :: distribution_function
    public :: distribution_property
    public :: distribution_recenter
    public :: multivariate_distribution
    public :: multivariate_distribution_function

    real(real64), parameter :: distribution_pi = 2.0d0 * acos(0.0d0)

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
    !! Compatibility fallback; override to avoid PDF underflow.
    class(distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst

    rst = probability_log(this%pdf(x))
end function

pure elemental function dist_survival(this, x) result(rst)
    !! Compatibility fallback; override to avoid cancellation in upper tails.
    class(distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst

    rst = 1.0d0 - this%cdf(x)
end function

pure elemental function dist_log_cdf(this, x) result(rst)
    !! Compatibility fallback; override for accurate extreme lower tails.
    class(distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst

    rst = probability_log(this%cdf(x))
end function

pure elemental function dist_log_survival(this, x) result(rst)
    !! Compatibility fallback; override for accurate extreme upper tails.
    class(distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst

    rst = probability_log(this%survival(x))
end function

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

end module fstats_distributions

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

module fstats_glm_types
    !! GLM family identifiers, solver controls, diagnostics, and vector callbacks.
    use iso_fortran_env
    implicit none
    private
    public :: FS_GLM_BINOMIAL
    public :: FS_GLM_POISSON
    public :: FS_GLM_GAMMA
    public :: glm_options
    public :: glm_convergence_info
    public :: glm_coefficient_statistics
    public :: inverse_link_function
    public :: robust_weight_function
    public :: assignment(=)

    integer(int32), parameter :: FS_GLM_BINOMIAL = 1
        !! Bernoulli/proportion responses with the logit link; no trial counts.
    integer(int32), parameter :: FS_GLM_POISSON = 2
        !! Nonnegative responses with the log link and V(mu) = mu.
    integer(int32), parameter :: FS_GLM_GAMMA = 3
        !! Positive responses with the reciprocal link and V(mu) = mu**2.

    type glm_options
        !! GLM controls. Component initialization and set_to_default agree.
        integer(int32) :: max_iteration_count = 100
            !! Positive maximum number of completed coefficient updates.
        real(real64) :: tolerance = 1.0d-8
            !! Finite positive absolute tolerance on maximum coefficient change.
            !! This is not a deviance or relative tolerance.
        logical :: use_robust_weighting = .true.
            !! Enable robust weights; set false for ordinary likelihood IRLS.
        real(real64) :: robust_weighting_constant = 4.685d0
            !! Finite positive cutoff when robust weighting is enabled.
            !! Residuals are raw y - mu, so this has response units. No scale
            !! estimate or residual standardization is performed.
    contains
        procedure, public :: set_to_default => go_set_defaults
            !! Restore all controls to their component defaults.
    end type

    interface assignment(=)
        !! Copy every component of a glm_options object.
        module procedure :: go_assign
    end interface

    type glm_coefficient_statistics
        !! Inferential statistics for one GLM coefficient, on the coefficient
        !! (linear-predictor) scale rather than the response scale. Use one
        !! object per fitted coefficient, in the same order as beta.
        !!
        !! This type is a container: components have no default values and
        !! must be populated by an inference calculation before use. IRLS can
        !! populate it through its optional stats output after convergence;
        !! on nonconvergence that output contains NaNs. Standard errors require a GLM covariance estimate;
        !! robust fits generally require an estimating-equation sandwich
        !! covariance, not simply the inverse robust-weighted information.
        real(real64) :: standard_error
            !! Square root of the coefficient covariance diagonal entry.
        real(real64) :: wald_statistic
            !! Signed coefficient divided by its standard error, testing the
            !! null hypothesis that the coefficient is zero. Ordinarily an
            !! asymptotic z statistic; approximate t inference may be used
            !! when dispersion is estimated. The reference distribution must
            !! be specified by the calculation producing this object.
        real(real64) :: p_value
            !! Two-sided tail probability under the null hypothesis and the
            !! chosen reference distribution, in [0,1]. This is not the
            !! probability that the coefficient is unimportant or zero.
        real(real64) :: confidence_interval_lower
            !! Lower confidence bound for the coefficient, not an interval
            !! half-width. The producing calculation specifies the confidence
            !! level and interval method (for example, Wald or profile likelihood).
        real(real64) :: confidence_interval_upper
            !! Upper confidence bound at the same level and using the same
            !! method as confidence_interval_lower. Separate bounds allow
            !! asymmetric intervals; no inverse-link transformation is implied.
    end type

    type glm_convergence_info
        !! Provides information regarding convergence status.
        logical :: converged = .false.
            !! True if the iteration process converged to the requested 
            !! tolerances; else, false.
        real(real64) :: residual_value = 0.0d0
            !! Maximum absolute coefficient change in the last update, not a
            !! response residual. IRLS initializes this to NaN before updates.
        integer(int32) :: iteration_count = 0
            !! Completed updates, excluding domain-preserving step halvings.
    end type

    abstract interface
        function inverse_link_function(eta) result(mu)
            !! Evaluate g**(-1)(eta) elementwise.
            use iso_fortran_env, only : real64
            real(real64), intent(in), dimension(:) :: eta
                !! N finite linear predictors in the inverse-link domain.
            real(real64), allocatable, dimension(:) :: mu
                !! N fitted means in the response family's mean domain.
        end function

        function robust_weight_function(r, c) result(w)
            !! Evaluate multiplicative robust weights on raw residuals.
            use iso_fortran_env, only : real64
            real(real64), intent(in), dimension(:) :: r
                !! N raw response residuals y - mu, not standardized.
            real(real64), intent(in) :: c
                !! Finite positive tuning cutoff in response units.
            real(real64), allocatable, dimension(:) :: w
                !! Allocated N-element vector of finite nonnegative weights.
                !! Zeros omit observations; remaining rows must identify beta.
        end function
    end interface

contains
! ------------------------------------------------------------------------------
    subroutine go_set_defaults(this)
        !! Sets the defaults for the [[glm_options]] object.
        class(glm_options), intent(inout) :: this
            !! The [[glm_options]] object.

        this%max_iteration_count = 100
        this%tolerance = 1.0d-8
        this%use_robust_weighting = .true.
        this%robust_weighting_constant = 4.685d0
    end subroutine

! ------------------------------------------------------------------------------
    subroutine go_assign(lhs, rhs)
        !! Assigns the contents of one [[glm_options]] object to another.
        type(glm_options), intent(inout) :: lhs
            !! The [[glm_options]] object to receive the contents of rhs.
        type(glm_options), intent(in) :: rhs
            !! The [[glm_options]] object to copy.

        lhs%max_iteration_count = rhs%max_iteration_count
        lhs%tolerance = rhs%tolerance
        lhs%use_robust_weighting = rhs%use_robust_weighting
        lhs%robust_weighting_constant = rhs%robust_weighting_constant
    end subroutine

! ------------------------------------------------------------------------------
end module
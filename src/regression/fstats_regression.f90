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

module fstats_regression
    use iso_fortran_env
    use ieee_arithmetic
    use linalg
    use fstats_errors
    use blas
    use fstats_descriptive_statistics
    use fstats_distributions
    use fstats_special_functions
    use fstats_hypothesis
    implicit none
    private
    public :: iteration_controls
    public :: convergence_info
    public :: regression_function
    public :: iteration_update
    public :: regression_statistics
    public :: r_squared
    public :: adjusted_r_squared
    public :: correlation
    public :: design_matrix
    public :: covariance_matrix
    public :: regression_covariance
    public :: calculate_regression_statistics
    public :: jacobian
    
! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    type regression_statistics
       !! A container for regression-related statistical information. 
        real(real64) :: standard_error
            !! The standard error for the model coefficient.
            !!
            !! $$ E_{s}(\beta_{i}) = \sqrt{\sigma^{2} C_{ii}} $$
        real(real64) :: t_statistic
            !! The T-statistic for the model coefficient.
            !!
            !! $$ t_o = \frac{ \beta_{i} }{E_{s}(\beta_{i})} $$
        real(real64) :: probability
            !! The probability that the coefficient is not statistically 
            !! important.  A statistically important coefficient will have a 
            !! low probability (p-value), typically 0.05 or lower; however, a 
            !! p-value of up to ~0.2 may be acceptable dependent upon the 
            !! problem.  Typically any p-value larger than ~0.2 indicates the 
            !! parameter is not statistically important for the model.
            !!
            !! $$ p = t_{|t_o|, df_{residual}} $$
        real(real64) :: confidence_interval
            !! The confidence interval for the parameter at the level 
            !! determined by the regression process.
            !!
            !! $$ c = t_{\alpha, df} E_{s}(\beta_{i}) $$
    end type

    type iteration_controls
        !! Provides a collection of iteration control parameters.
        integer(int32) :: max_iteration_count
            !! Defines the maximum number of iterations allowed.
        integer(int32) :: max_function_evaluations
            !! Defines the maximum number of function evaluations allowed.
        real(real64) :: gradient_tolerance
            !! Defines a tolerance on the gradient of the fitted function.
        real(real64) :: change_in_solution_tolerance
            !! Defines a tolerance on the change in parameter values.
        real(real64) :: residual_tolerance
            !! Defines a tolerance on the metric associated with the residual 
            !! error.
        real(real64) :: iteration_improvement_tolerance
            !! Defines a tolerance to ensure adequate improvement on each 
            !! iteration.
        integer(int32) :: max_iteration_between_updates
            !! Defines how many iterations can pass before a re-evaluation of 
            !! the Jacobian matrix is forced.
    contains
        procedure, public :: set_to_default => lm_set_default_tolerances
    end type

    type convergence_info
        !! Provides information regarding convergence status.
        logical :: converge_on_gradient
            !! True if convergence on the gradient was achieved; else, false.
        real(real64) :: gradient_value
            !! The value of the gradient test parameter.
        logical :: converge_on_solution_change
            !! True if convergence on the change in solution was achieved; else,
            !! false.
        real(real64) :: solution_change_value
            !! The value of the change in solution parameter.
        logical :: converge_on_residual_parameter
            !! True if convergence on the residual error parameter was achieved; 
            !! else, false.
        real(real64) :: residual_value
            !! The value of the residual error parameter.
        logical :: reach_iteration_limit
            !! True if the solution did not converge in the allowed number of 
            !! iterations.
        integer(int32) :: iteration_count
            !! The iteration count.
        logical :: reach_function_evaluation_limit
            !! True if the solution did not converge in the allowed number of
            !! function evaluations.
        integer(int32) :: function_evaluation_count
            !! The function evaluation count.
        logical :: user_requested_stop
            !! True if the user requested the stop; else, false.
        integer(int32) :: covariance_rank = 0
            !! Numerical rank of the final weighted Jacobian when covariance is requested.
        real(real64) :: covariance_condition_number = 0.0d0
            !! Condition number of that Jacobian; infinity indicates rank deficiency.
    end type

    interface
        subroutine regression_function(xdata, params, f, stop, args)
            !! Defines the interface of a subroutine computing the function
            !! values at each of the N data points as part of a regression
            !! analysis.
            use iso_fortran_env, only : real64
            real(real64), intent(in), dimension(:) :: xdata
                !! An N-element array containing the N independent data points.
            real(real64), intent(in), dimension(:) :: params
                !! An M-element array containing the M model parameters.
            real(real64), intent(out), dimension(:) :: f
                !! An N-element array where the results of the N function 
                !! evaluations will be written.
            logical, intent(out) :: stop
                !! A mechanism to force a stop to the iteration process.  If
                !! set to true, the iteration process will terminate.  If set
                !! to false, the iteration process will continue along as 
                !! normal.
            class(*), intent(inout), optional :: args
                !! An optional argument allowing the passing in/out of data.
        end subroutine

        subroutine iteration_update(iter, funvals, resid, params, step)
            !! Defines a routine for providing updates about an iteration
            !! process.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: iter
                !! The current iteration number.
            real(real64), intent(in), dimension(:) :: funvals
                !! The function values.
            real(real64), intent(in), dimension(:) :: resid
                !! The residuals.
            real(real64), intent(in), dimension(:) :: params
                !! The model parameters.
            real(real64), intent(in), dimension(:) :: step
                !! Step sizes for each parameter.
        end subroutine
    end interface

contains

! ------------------------------------------------------------------------------
pure function r_squared(x, xm) result(rst)
    !! Computes the R-squared value for a data set.
    !!
    !! The R-squared value is computed by determining the sum of the squares
    !! of the residuals: 
    !! $$ SS_{res} = \Sigma \left( y_i - f_i \right)^2 $$
    !! The total sum of the squares: 
    !! $$ SS_{tot} = \Sigma \left( y_i - \bar{y} \right)^2 $$. 
    !! The R-squared value is then: 
    !! $$ R^2 = 1 - \frac{SS_{res}}{SS_{tot}} $$.
    !!
    !! See Also:
    !!
    !! - <a href="https://en.wikipedia.org/wiki/Coefficient_of_determination" target="_blank">Wikipedia</a>
    real(real64), intent(in) :: x(:)
        !! An N-element array containing the dependent variables from 
        !! the data set.
    real(real64), intent(in) :: xm(:)
        !! An N-element array containing the corresponding modeled 
        !! values.
    real(real64) :: rst
        !! The result.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: esum, vt
    
    ! Initialization

    ! Input Check
    n = size(x)
    if (size(xm) /= n) error stop FS_ARRAY_SIZE_ERROR

    ! Process
    esum = zero
    do i = 1, n
        esum = esum + (x(i) - xm(i))**2
    end do
    vt = variance(x) * (n - one)
    rst = one - esum / vt
end function

! ------------------------------------------------------------------------------
pure function adjusted_r_squared(p, x, xm) result(rst)
    !! Computes the adjusted R-squared value for a data set.
    !!
    !! The adjusted R-squared provides a mechanism for tempering the effects
    !! of extra explanatory variables on the traditional R-squared 
    !! calculation.  It is computed by noting the sample size \( n \) and 
    !! the number of variables \( p \).
    !! $$ \bar{R}^2 = 1 - \left( 1 - R^2 \right) \frac{n - 1}{n - p} $$.
    !!
    !! See Also:
    !!
    !! - <a href="https://en.wikipedia.org/wiki/Coefficient_of_determination#Adjusted_R2" target="_blank">Wikipedia</a>
    integer(int32), intent(in) :: p
        !! The number of variables.
    real(real64), intent(in) :: x(:)
        !! An N-element array containing the dependent variables from 
        !! the data set.
    real(real64), intent(in) :: xm(:)
        !! An N-element array containing the corresponding modeled 
        !! values.
    real(real64) :: rst
        !! The result.

    ! Local Variables
    integer(int32) :: n
    real(real64) :: r2

    ! Parameters
    real(real64), parameter :: one = 1.0d0
    
    ! Initialization
    n = size(x)

    ! Process
    r2 = r_squared(x, xm)
    rst = one - (one - r2) * (n - one) / (n - p - one)
end function

! ------------------------------------------------------------------------------
pure function correlation(x, y) result(rst)
    !! Computes the sample correlation coefficient (an estimate to the 
    !! population Pearson correlation) as follows.
    !!
    !! $$ r_{xy} = \frac{cov(x, y)}{s_{x} s_{y}} $$.
    !!
    !! Where, \( s_{x} \) & \( s_{y} \) are the sample standard deviations of
    !! x and y respectively.
    real(real64), intent(in), dimension(:) :: x
        !! The first N-element data set.
    real(real64), intent(in), dimension(size(x)) :: y
        !! The second N-element data set.
    real(real64) :: rst
        !! Correlation coefficient on [-1, 1]; constant, nonfinite, or
        !! undersized data give NaN. Scaled centering avoids overflowing moments.

    real(real64), allocatable, dimension(:) :: centeredX, centeredY
    real(real64) :: scaleX, scaleY, normX, normY
    integer(int32) :: observation

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (size(x) < 2 .or. .not.all(ieee_is_finite(x)) .or. .not.all(ieee_is_finite(y))) return
    scaleX = maxval(abs(x))
    scaleY = maxval(abs(y))
    if (scaleX == 0.0d0 .or. scaleY == 0.0d0) return
    allocate(centeredX(size(x)), centeredY(size(x)))
    do observation = 1, size(x)
        if ((x(observation) >= 0.0d0 .and. x(1) >= 0.0d0) .or. &
            (x(observation) <= 0.0d0 .and. x(1) <= 0.0d0)) then
            centeredX(observation) = (x(observation) - x(1)) / scaleX
        else
            centeredX(observation) = x(observation) / scaleX - x(1) / scaleX
        end if
        if ((y(observation) >= 0.0d0 .and. y(1) >= 0.0d0) .or. &
            (y(observation) <= 0.0d0 .and. y(1) <= 0.0d0)) then
            centeredY(observation) = (y(observation) - y(1)) / scaleY
        else
            centeredY(observation) = y(observation) / scaleY - y(1) / scaleY
        end if
    end do
    centeredX = centeredX - mean(centeredX)
    centeredY = centeredY - mean(centeredY)
    normX = norm2(centeredX)
    normY = norm2(centeredY)
    if (normX == 0.0d0 .or. normY == 0.0d0) return
    rst = max(-1.0d0, min(1.0d0, dot_product(centeredX / normX, centeredY / normY)))
end function

! ------------------------------------------------------------------------------
pure function design_matrix(order, intercept, x) result(c)
    !! Computes the design matrix \( X \) for the linear 
    !! least-squares regression problem of \( X \beta = y \), where 
    !! \( X \) is the matrix computed here, \( \beta \) is 
    !! the vector of coefficients to be determined, and \( y \) is the 
    !! vector of measured dependent variables.
    !!
    !! See Also
    !!
    !! - <a href="https://en.wikipedia.org/wiki/Linear_regression" target="_blank">Wikipedia - Linear Regression</a>
    !! - <a href="https://en.wikipedia.org/wiki/Vandermonde_matrix" target="_blank">Wikipedia - Vandermonde Matrix</a>
    !! - <a href="https://en.wikipedia.org/wiki/Design_matrix" target="_blank">Wikipedia - Design Matrix</a>
    integer(int32), intent(in) :: order
        !! The order of the equation to fit.  This value must be
        !! at least one (linear equation), but can be higher as desired.
    logical, intent(in) :: intercept
        !! Set to true if the intercept is being computed
        !! as part of the regression; else, false.
    real(real64), intent(in) :: x(:)
        !! An N-element array containing the independent variable
        !! measurement points.
    real(real64), allocatable :: c(:,:)
        !! An N-by-K matrix where the results will be written.  K
        !! must equal order + 1 in the event intercept is true; 
        !! however, if intercept is false, K must equal order.

    ! Parameters
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, start, npts, ncols
    
    ! Initialization
    npts = size(x)
    ncols = order
    if (intercept) ncols = ncols + 1
    allocate(c(npts, ncols))

    ! Input Check
    if (order < 1) error stop FS_INVALID_INPUT_ERROR

    ! Process
    if (intercept) then
        c(:,1) = one
        c(:,2) = x
        start = 3
    else
        c(:,1) = x
        start = 2
    end if
    if (start > ncols) return
    do i = start, ncols
        c(:,i) = c(:,i-1) * x
    end do
end function

! ------------------------------------------------------------------------------
pure function covariance_matrix(x) result(c)
    !! Computes the covariance matrix \( C \) where 
    !! \( C = \left( X^{T} X \right)^{-1} \) and \( X \) is computed
    !! by design_matrix.
    !! Uses an SVD of X, without forming X**T * X. Rank-deficient inputs give
    !! the Moore-Penrose covariance; use regression_covariance for rank diagnostics.
    !!
    !! See Also
    !!
    !! - <a href="https://en.wikipedia.org/wiki/Covariance_matrix" target="_blank">Wikipedia - Covariance Matrix</a>
    !! - <a href="https://en.wikipedia.org/wiki/Linear_regression" target="_blank">Wikipedia - Linear Regression</a>
    real(real64), intent(in) :: x(:,:)
        !! An M-by-N matrix containing the formatted independent data
        !!  matrix \( X \) as computed by design_matrix.
    real(real64), allocatable :: c(:,:)
        !! The N-by-N covariance matrix.

    call regression_covariance(x, c)
end function

! ------------------------------------------------------------------------------
pure subroutine regression_covariance(x, c, rank, condition_number)
    !! Computes the unscaled regression covariance using an SVD of X itself.
    !! Singular directions below MAX(M, N) * EPSILON * largest singular value
    !! are discarded, giving the Moore-Penrose covariance on identifiable directions.
    !! Rank-deficient fits require care: a zero pseudocovariance along a discarded
    !! direction does not imply that the corresponding parameter is known exactly.
    real(real64), intent(in), dimension(:,:) :: x
        !! M-by-N design matrix, or a Jacobian whose rows include sqrt(weights).
    real(real64), intent(out), allocatable, dimension(:,:) :: c
        !! N-by-N covariance before residual-variance scaling. Nonfinite X gives NaNs.
    integer(int32), intent(out), optional :: rank
        !! Number of retained singular directions; zero for empty/invalid input.
    real(real64), intent(out), optional :: condition_number
        !! Condition number of X; +infinity for rank deficiency, NaN for invalid X.
    real(real64), allocatable, dimension(:) :: singular_values
    real(real64), allocatable, dimension(:,:) :: vectors, scaled_vectors
    real(real64) :: relative_tolerance
    integer(int32) :: ncols, direction, effective_rank

    ncols = size(x, 2)
    allocate(c(ncols, ncols), source = 0.0d0)
    if (present(rank)) rank = 0
    if (present(condition_number)) condition_number = ieee_value(0.0d0, ieee_positive_inf)
    if (.not.all(ieee_is_finite(x))) then
        c = ieee_value(0.0d0, ieee_quiet_nan)
        if (present(condition_number)) condition_number = ieee_value(0.0d0, ieee_quiet_nan)
        return
    end if
    if (minval(shape(x)) == 0) return
    call svd(x, s = singular_values, vt = vectors)
    if (singular_values(1) == 0.0d0) return
    relative_tolerance = real(maxval(shape(x)), real64) * epsilon(1.0d0)
    allocate(scaled_vectors(size(singular_values), ncols), source = 0.0d0)
    effective_rank = 0
    do direction = 1, size(singular_values)
        if (singular_values(direction) / singular_values(1) <= relative_tolerance) cycle
        effective_rank = effective_rank + 1
        scaled_vectors(direction,:) = vectors(direction,:) / singular_values(direction)
    end do
    c = matmul(transpose(scaled_vectors), scaled_vectors)
    if (present(rank)) rank = effective_rank
    if (present(condition_number)) then
        if (effective_rank == ncols) &
            condition_number = singular_values(1) / singular_values(ncols)
    end if
end subroutine

! ------------------------------------------------------------------------------
function calculate_regression_statistics(resid, params, c, alpha) &
    result(rst)
    !! Computes statistics for the quality of fit for a regression 
    !! model.
    real(real64), intent(in) :: resid(:)
        !! An M-element array containing the model residual errors.
    real(real64), intent(in) :: params(:)
        !! An N-element array containing the model parameters.
    real(real64), intent(in) :: c(:,:)
        !! The N-by-N covariance matrix.
    real(real64), intent(in), optional :: alpha
        !! The significance level at which to evaluate the confidence 
        !! intervals.  The default value is 0.05 such that a 95% 
        !! confidence interval is calculated.
    type(regression_statistics), allocatable :: rst(:)
        !! A regression_statistics object containing the analysis results.

    ! Parameters
    real(real64), parameter :: p05 = 0.05d0
    real(real64), parameter :: half = 0.5d0
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, m, n, dof
    real(real64) :: a, ssr, var, talpha
    type(t_distribution) :: dist

    ! Initialization
    m = size(resid)
    n = size(params)
    dof = m - n
    if (present(alpha)) then
        a = alpha
    else
        a = p05
    end if

    ! Input Checking
    if (dof <= 0) error stop FS_INVALID_INPUT_ERROR
    if (a <= zero .or. a >= one) error stop FS_INVALID_INPUT_ERROR
    if (size(c, 1) /= n .or. size(c, 2) /= n) error stop FS_MATRIX_SIZE_ERROR
    allocate(rst(n))

    ! Process
    ssr = norm2(resid)**2   ! sum of the squares of the residual
    var = ssr / dof
    dist%dof = real(dof, real64)
    talpha = confidence_interval(dist, a, one, 1)
    do i = 1, n
        rst(i)%standard_error = sqrt(var * c(i,i))
        rst(i)%confidence_interval = talpha * rst(i)%standard_error
        if (ssr == zero) then
            if (params(i) == zero) then
                rst(i)%t_statistic = zero
                rst(i)%probability = one
            else
                rst(i)%t_statistic = sign( &
                    ieee_value(params(i), IEEE_POSITIVE_INF), params(i) &
                )
                rst(i)%probability = zero
            end if
        else
            rst(i)%t_statistic = params(i) / rst(i)%standard_error
            rst(i)%probability = regularized_beta( &
                half * dof, &
                half, &
                real(dof, real64) / (dof + (rst(i)%t_statistic)**2) &
            )
        end if
    end do
end function

! ------------------------------------------------------------------------------
subroutine jacobian(fun, xdata, params, &
    jac, stop, f0, f1, step, args)
    !! Computes the Jacobian matrix for a nonlinear regression problem.
    procedure(regression_function), intent(in), pointer :: fun
        !! A pointer to the regression_function to evaluate.
    real(real64), intent(in) :: xdata(:)
        !! The M-element array containing x-coordinate data.
    real(real64), intent(in) :: params(:)
        !! The N-element array containing the model parameters.
    real(real64), intent(out) :: jac(:,:)
        !! The M-by-N matrix where the Jacobian will be written.
    logical, intent(out) :: stop
        !! A value that the user can set in fun forcing the
        !! evaluation process to stop prior to completion.
    real(real64), intent(in), optional, target :: f0(:)
        !! An optional M-element array containing the model values
        !!  using the current parameters as defined in m.  This input 
        !! can be used to prevent the routine from performing a 
        !! function evaluation at the model parameter state defined in 
        !! params.
    real(real64), intent(out), optional, target :: f1(:)
        !! An optional M-element workspace array used for function
        !! evaluations.
    real(real64), intent(in), optional :: step
        !! The differentiation step size.  The default is the square 
        !! root of machine precision.
    class(*), intent(inout), optional :: args
        !! An optional argument allowing the passing in/out of data.

    ! Local Variables
    real(real64) :: h
    integer(int32) :: m, n, expected, actual
    real(real64), pointer :: f1p(:), f0p(:)
    real(real64), allocatable, target :: f1a(:), f0a(:), work(:)

    ! Initialization
    if (present(step)) then
        h = step
    else
        h = sqrt(epsilon(h))
    end if
    m = size(xdata)
    n = size(params)

    ! Input Size Checking
    if (size(jac, 1) /= m .or. size(jac, 2) /= n) error stop FS_MATRIX_SIZE_ERROR
    if (present(f0)) then
        ! Check Size
        if (size(f0) /= m) error stop FS_ARRAY_SIZE_ERROR
        f0p(1:m) => f0
    else
        ! Allocate space, and fill the array with the current function
        ! results
        allocate(f0a(m))
        f0p(1:m) => f0a
        call fun(xdata, params, f0p, stop, args = args)
        if (stop) return
    end if
    if (present(f1)) then
        ! Check Size
        if (size(f1) /= m) error stop FS_ARRAY_SIZE_ERROR
        f1p(1:m) => f1
    else
        ! Allocate space
        allocate(f1a(m))
        f1p(1:m) => f1a
    end if

    ! Allocate a workspace array the same size as params
    allocate(work(n))

    ! Compute the Jacobian
    call jacobian_finite_diff(fun, xdata, params, f0p, jac, f1p, &
        stop, h, work, args = args)
end subroutine

! ******************************************************************************
! PRIVATE ROUTINES
! ------------------------------------------------------------------------------
! Computes the Jacobian matrix via a forward difference.
!
! Inputs:
! - fun: The function to evaluate
! - xdata: The independent coordinate data to fit (M-by-1)
! - params: The model parameters (N-by-1)
! - f0: The current model estimate (M-by-1)
! - step: The differentiation step size
!
! Outputs:
! - jac: The Jacobian matrix (M-by-N)
! - f1: A workspace array for the model output (M-by-1)
! - stop: A flag allowing the user to terminate model execution
! - work: A workspace array for the model parameters (N-by-1)
subroutine jacobian_finite_diff(fun, xdata, params, f0, jac, f1, &
    stop, step, work, args)
    ! Arguments
    procedure(regression_function), intent(in), pointer :: fun
    real(real64), intent(in), dimension(:) :: xdata, params
    real(real64), intent(in), dimension(:) :: f0
    real(real64), intent(out), dimension(:,:) :: jac
    real(real64), intent(out), dimension(:) :: f1, work
    logical, intent(out) :: stop
    real(real64), intent(in) :: step
    class(*), intent(inout), optional :: args

    ! Local Variables
    integer(int32) :: i, n

    ! Initialization
    n = size(params)

    ! Cycle over each column of the Jacobian and calculate the derivative
    ! via a forward difference scheme
    !
    ! J(i,j) = df(i) / dx(j)
    work = params
    do i = 1, n
        work(i) = work(i) + step
        call fun(xdata, work, f1, stop, args = args)
        if (stop) return

        jac(:,i) = (f1 - f0) / step
        work(i) = params(i)
    end do
end subroutine

! ******************************************************************************
! SETTINGS DEFAULTS
! ------------------------------------------------------------------------------
! Sets up default tolerances.
subroutine lm_set_default_tolerances(x)
    ! Arguments
    class(iteration_controls), intent(inout) :: x

    ! Set defaults
    x%max_iteration_count = 500
    x%max_function_evaluations = 5000
    x%max_iteration_between_updates = 10
    x%gradient_tolerance = 1.0d-8
    x%residual_tolerance = 0.5d-2
    x%change_in_solution_tolerance = 1.0d-6
    x%iteration_improvement_tolerance = 1.0d-1
end subroutine

! ------------------------------------------------------------------------------
end module
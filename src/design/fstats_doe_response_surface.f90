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

module fstats_doe_response_surface
    !! Response-surface optimization.
    !!
    !! This module locates the factor settings that maximize, minimize, or 
    !! achieve a target value of a fitted doe_model response surface within a
    !! box of factor bounds.  The search is performed in coded variables on 
    !! the range [-1, 1], which is the scale on which the model is assumed 
    !! to have been fit, using a projected gradient method with a 
    !! backtracking line search that is started from multiple points to 
    !! reduce the chance of stopping at a local optimum.
    use iso_fortran_env, only : int32, real64
    use fstats_errors, only : FS_INVALID_INPUT_ERROR, FS_ARRAY_SIZE_ERROR
    use fstats_experimental_design, only : doe_model, doe_optimization_result
    use fstats_doe_models, only : doe_evaluate_model
    use linalg, only : eigen
    implicit none
    private

    public :: doe_optimize_rsm
    public :: doe_stationary_point
contains
! ------------------------------------------------------------------------------
subroutine doe_stationary_point(mdl, n, z, nature, ok, eigenvalues)
    !! Computes the stationary point of a quadratic (nway <= 2) response 
    !! surface in coded variables and classifies its nature.
    !!
    !! The model is written as
    !!
    !! $$ \hat{y}(z) = \beta_{0} + g^{T} z + \frac{1}{2} z^{T} H z $$
    !!
    !! where the gradient \(g\) and Hessian \(H\) are recovered exactly from
    !! model evaluations.  The stationary point solves \(H z = -g\), and is 
    !! classified by the eigenvalues of \(H\): all negative is a maximum, 
    !! all positive is a minimum, mixed signs is a saddle, and a (near) zero 
    !! eigenvalue is singular (a ridge or flat direction, with no unique 
    !! stationary point).
    class(doe_model), intent(in) :: mdl
        !! The fitted DOE model, defined in coded variables.  Its nway 
        !! must be 1 or 2.
    integer, intent(in) :: n
        !! The number of factors.
    real(real64), allocatable, intent(out), dimension(:) :: z
        !! The stationary point in coded variables.  This is allocated only 
        !! if ok is true.
    character(len=:), allocatable, intent(out) :: nature
        !! One of `maximum`, `minimum`, `saddle`, `singular`, or 
        !! `not quadratic`.
    logical, intent(out) :: ok
        !! True if a unique stationary point was found.
    real(real64), allocatable, intent(out), dimension(:), optional :: &
        eigenvalues
        !! The eigenvalues of the Hessian, in coded variables.

    integer :: i, j
    real(real64) :: y0, lmax
    real(real64), allocatable, dimension(:) :: g, ev, ep, em, e0
    real(real64), allocatable, dimension(:,:) :: h, vecs
    
    ok = .false.
    if (mdl%nway < 1 .or. mdl%nway > 2 .or. n < 1) then
        nature = "not quadratic"
        return
    end if

    ! Recover g and H exactly (the model is quadratic)
    allocate(g(n), ev(n), ep(n), em(n), e0(n), h(n, n))
    e0 = 0.0d0
    y0 = response(mdl, e0)
    do i = 1, n
        e0 = 0.0d0
        e0(i) = 1.0d0
        ep(i) = response(mdl, e0)
        e0(i) = -1.0d0
        em(i) = response(mdl, e0)
        g(i) = 0.5d0 * (ep(i) - em(i))
        h(i, i) = ep(i) + em(i) - 2.0d0 * y0
    end do
    do i = 1, n - 1
        do j = i + 1, n
            e0 = 0.0d0
            e0(i) = 1.0d0
            e0(j) = 1.0d0
            h(i, j) = response(mdl, e0) - ep(i) - ep(j) + y0
            h(j, i) = h(i, j)
        end do
    end do

    ! Classify
    call eigen(h, ev, vecs)
    if (present(eigenvalues)) then
        allocate(eigenvalues, source = ev)
    end if
    lmax = maxval(abs(ev))
    if (lmax == 0.0d0 .or. minval(abs(ev)) <= 1.0d-10 * lmax) then
        nature = "singular"
        return
    end if
    if (all(ev < 0.0d0)) then
        nature = "maximum"
    else if (all(ev > 0.0d0)) then
        nature = "minimum"
    else
        nature = "saddle"
    end if

    ! Solve H z = -g using the spectral decomposition H = V L V**T
    allocate(z(n))
    z = -matmul(vecs, matmul(transpose(vecs), g) / ev)
    ok = .true.
end subroutine

! ------------------------------------------------------------------------------
function doe_optimize_rsm(mdl, x_low, x_high, goal, target, tol, max_iter) &
    result(opt)
    !! Optimizes a fitted response surface within a box of factor bounds.
    !!
    !! The model is assumed to be defined in coded variables 
    !! \(z_{i} \in [-1, 1]\), which are related to the natural variables by
    !!
    !! $$ x_{i} = \frac{x_{high,i} + x_{low,i}}{2} + 
    !! \frac{x_{high,i} - x_{low,i}}{2} z_{i} $$
    !!
    !! The problem solved depends on the goal:
    !!
    !! - `max`: maximize \(\hat{y}(z)\)
    !! - `min`: minimize \(\hat{y}(z)\)
    !! - `target`: minimize \((\hat{y}(z) - y_{target})^{2}\)
    !!
    !! subject to \(-1 \le z_{i} \le 1\).  A projected gradient method with 
    !! an Armijo backtracking line search is used.  Gradients are computed 
    !! using central differences, which are essentially exact for the 
    !! polynomial models used here.  The search is repeated from the design 
    !! center and from points near the corners of the coded region (up to 
    !! 64 corner patterns), and the best result is returned.  A run is 
    !! considered converged when the norm of the projected gradient step, 
    !! \(\| z - P(z - \nabla f) \|\), falls below the tolerance.
    !!
    !! For quadratic models (nway <= 2) and the `max` or `min` goals, the 
    !! stationary point is first computed analytically (see 
    !! doe_stationary_point).  If it lies within the bounds and is of the 
    !! matching type (a maximum for `max`, a minimum for `min`), it is the 
    !! global optimum and is returned directly.  Otherwise the iterative 
    !! search is used.
    !!
    !! Note, the result is a local optimum of the fitted model, found from 
    !! multiple starting points; it is not a guarantee of a global optimum.
    !! The optimum may lie on the boundary of the factor region, and 
    !! extrapolation beyond the region of the experimental data should be
    !! confirmed experimentally.
    class(doe_model), intent(in) :: mdl
        !! The fitted DOE model, defined in coded variables.
    real(real64), intent(in), dimension(:) :: x_low
        !! The lower natural-scale bound for each factor.
    real(real64), intent(in), dimension(:) :: x_high
        !! The upper natural-scale bound for each factor.  This array must be
        !! the same size as x_low, and each element must exceed the 
        !! corresponding element of x_low.
    character(len=*), intent(in), optional :: goal
        !! The optimization goal: `max` (the default), `min`, or `target`.
    real(real64), intent(in), optional :: target
        !! The target response value.  This is required if goal is `target`,
        !! and is ignored otherwise.
    real(real64), intent(in), optional :: tol
        !! The convergence tolerance on the projected gradient step.  The 
        !! default is 1.0d-8.
    integer, intent(in), optional :: max_iter
        !! The maximum number of iterations allowed for each starting point.
        !! The default is 1000.
    type(doe_optimization_result) :: opt
        !! A doe_optimization_result containing the optimal factors in both 
        !! coded and natural scale, the predicted response at the optimum, 
        !! the iteration count and convergence flag of the best run, the 
        !! method name, and the tolerance used.

    ! Parameters
    integer, parameter :: max_corner_bits = 6

    ! Local Variables
    integer :: i, k, n, nbits, nstarts, mi, iters, best_iters, kind
    logical :: ok, best_ok, have_best
    real(real64) :: tl, tgt, f, best_f
    real(real64), allocatable, dimension(:) :: z, z0, zbest, zs
    character(len=:), allocatable :: nature

    ! Initialization
    n = size(x_low)
    if (size(x_high) /= n) error stop FS_ARRAY_SIZE_ERROR
    if (n < 1) error stop FS_INVALID_INPUT_ERROR
    if (any(x_high <= x_low)) error stop FS_INVALID_INPUT_ERROR
    tl = 1.0d-8
    if (present(tol)) tl = tol
    if (tl <= 0.0d0) error stop FS_INVALID_INPUT_ERROR
    mi = 1000
    if (present(max_iter)) mi = max_iter
    if (mi < 1) error stop FS_INVALID_INPUT_ERROR

    kind = 1
    if (present(goal)) then
        select case (trim(adjustl(goal)))
        case ("max")
            kind = 1
        case ("min")
            kind = 2
        case ("target")
            kind = 3
        case default
            error stop FS_INVALID_INPUT_ERROR
        end select
    end if
    tgt = 0.0d0
    if (kind == 3) then
        if (.not.present(target)) error stop FS_INVALID_INPUT_ERROR
        tgt = target
    end if

    allocate(z(n), z0(n), zbest(n))

    ! Analytic solution: an interior stationary point of the matching type 
    ! of a quadratic surface is the global optimum
    if (kind /= 3 .and. mdl%nway <= 2) then
        call doe_stationary_point(mdl, n, zs, nature, ok)
        if (ok) then
            if (all(abs(zs) <= 1.0d0) .and. &
                ((kind == 1 .and. nature == "maximum") .or. &
                (kind == 2 .and. nature == "minimum"))) then
                allocate(opt%optimal_coded_factors, source = zs)
                allocate(opt%optimal_natural_factors(n))
                opt%optimal_natural_factors = 0.5d0 * (x_high + x_low) &
                    + 0.5d0 * (x_high - x_low) * zs
                opt%optimal_response = response(mdl, zs)
                opt%converged = .true.
                opt%method = "analytic stationary point"
                opt%convergence_tolerance = tl
                opt%iteration_count = 0
                return
            end if
        end if
    end if

    ! Multi-start search: the center followed by corner patterns
    nbits = min(n, max_corner_bits)
    nstarts = 1 + 2**nbits
    have_best = .false.
    best_f = 0.0d0
    best_iters = 0
    best_ok = .false.
    do k = 0, nstarts - 1
        if (k == 0) then
            z0 = 0.0d0
        else
            do i = 1, n
                if (btest(k - 1, mod(i - 1, nbits))) then
                    z0(i) = 0.9d0
                else
                    z0(i) = -0.9d0
                end if
            end do
        end if
        call projected_gradient(mdl, kind, tgt, z0, tl, mi, z, f, iters, ok)
        if (.not.have_best .or. f < best_f) then
            have_best = .true.
            best_f = f
            zbest = z
            best_iters = iters
            best_ok = ok
        end if
    end do

    ! Output
    allocate(opt%optimal_coded_factors, source = zbest)
    allocate(opt%optimal_natural_factors(n))
    opt%optimal_natural_factors = 0.5d0 * (x_high + x_low) &
        + 0.5d0 * (x_high - x_low) * zbest
    opt%optimal_response = response(mdl, zbest)
    opt%converged = best_ok
    opt%method = "projected gradient (multi-start)"
    opt%convergence_tolerance = tl
    opt%iteration_count = best_iters
end function

! ------------------------------------------------------------------------------
function response(mdl, z) result(rst)
    !! Evaluates the model at a single coded point.
    class(doe_model), intent(in) :: mdl
    real(real64), intent(in), dimension(:) :: z
    real(real64) :: rst

    real(real64), allocatable, dimension(:) :: y
    
    y = doe_evaluate_model(mdl, reshape(z, [1, size(z)]))
    rst = y(1)
end function

! ------------------------------------------------------------------------------
function objective(mdl, kind, target, z) result(rst)
    !! Evaluates the objective (to be minimized) at a coded point.
    class(doe_model), intent(in) :: mdl
    integer, intent(in) :: kind
    real(real64), intent(in) :: target
    real(real64), intent(in), dimension(:) :: z
    real(real64) :: rst

    rst = response(mdl, z)
    select case (kind)
    case (1)
        rst = -rst
    case (3)
        rst = (rst - target)**2
    end select
end function

! ------------------------------------------------------------------------------
subroutine projected_gradient(mdl, kind, target, z0, tol, max_iter, z, f, &
    iters, converged)
    !! Minimizes the objective over the box [-1, 1]**N using a projected 
    !! gradient method with Armijo backtracking.
    class(doe_model), intent(in) :: mdl
    integer, intent(in) :: kind
    real(real64), intent(in) :: target
    real(real64), intent(in), dimension(:) :: z0
    real(real64), intent(in) :: tol
    integer, intent(in) :: max_iter
    real(real64), intent(out), dimension(:) :: z
    real(real64), intent(out) :: f
    integer, intent(out) :: iters
    logical, intent(out) :: converged

    ! Parameters
    real(real64), parameter :: c1 = 1.0d-4
    real(real64), parameter :: h = 1.0d-6
    real(real64), parameter :: min_step = 1.0d-16

    ! Local Variables
    integer :: i, n
    logical :: accepted
    real(real64) :: step, fnew, decrease
    real(real64), allocatable, dimension(:) :: g, znew

    n = size(z0)
    allocate(g(n), znew(n))
    z = z0
    f = objective(mdl, kind, target, z)
    step = 1.0d0
    converged = .false.
    iters = 0
    do while (iters < max_iter)
        call gradient(mdl, kind, target, z, h, g)

        ! Projected gradient step norm
        znew = max(-1.0d0, min(1.0d0, z - g))
        if (norm2(z - znew) <= tol) then
            converged = .true.
            exit
        end if

        ! Backtracking line search along the projected path
        accepted = .false.
        step = min(step * 2.0d0, 1.0d6)
        do while (step > min_step)
            znew = max(-1.0d0, min(1.0d0, z - step * g))
            fnew = objective(mdl, kind, target, znew)
            decrease = c1 * dot_product(g, z - znew)
            if (fnew <= f - decrease) then
                accepted = .true.
                exit
            end if
            step = 0.5d0 * step
        end do
        iters = iters + 1
        if (.not.accepted) then
            ! No further progress is possible at machine precision
            converged = norm2(z - max(-1.0d0, min(1.0d0, z - g))) <= &
                sqrt(tol)
            exit
        end if
        z = znew
        f = fnew
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine gradient(mdl, kind, target, z, h, g)
    !! Computes the objective gradient using central differences.  Points 
    !! are allowed to step slightly outside the bounds as the model is a 
    !! polynomial that is defined everywhere.
    class(doe_model), intent(in) :: mdl
    integer, intent(in) :: kind
    real(real64), intent(in) :: target
    real(real64), intent(in), dimension(:) :: z
    real(real64), intent(in) :: h
    real(real64), intent(out), dimension(:) :: g

    integer :: i
    real(real64), allocatable, dimension(:) :: zp, zm

    allocate(zp(size(z)), zm(size(z)))
    do i = 1, size(z)
        zp = z
        zm = z
        zp(i) = z(i) + h
        zm(i) = z(i) - h
        g(i) = (objective(mdl, kind, target, zp) &
            - objective(mdl, kind, target, zm)) / (2.0d0 * h)
    end do
end subroutine

end module fstats_doe_response_surface

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

module fstats_doe_models
    !! DOE model fitting and evaluation.
    use iso_fortran_env, only : int32, real64
    use fstats_errors, only : FS_INVALID_INPUT_ERROR, FS_ARRAY_SIZE_ERROR
    use fstats_regression, only : regression_statistics, covariance_matrix, &
        calculate_regression_statistics
    use fstats_experimental_design, only : doe_model
    implicit none
    private

    public :: doe_fit_model
    public :: doe_evaluate_model
    public :: doe_model_covariance
    public :: doe_design_terms

    interface doe_evaluate_model
        module procedure :: doe_evaluate_model_1
        module procedure :: doe_evaluate_model_2
    end interface
contains
! ------------------------------------------------------------------------------
function doe_fit_model(nway, x, y, map, alpha) result(rst)
    use blas, only : DGEMM
    use ieee_arithmetic
    !! Fits a Taylor series model to the provided data.
    !!
    !! $$ Y = \beta_{0} + \sum_{i=1}^{n} \beta_{i} X_{i} + \sum_{i=1}^{n-1} 
    !! \sum_{j=i+1}^{n} \beta_{ij} X_{i} X_{j} + \sum_{i=1}^{n} \beta_{ii} 
    !! X_{i}^{2} + \sum_{i=1}^{n-2} \sum_{j=i+1}^{n-1} \sum_{k=j+1}^{n} 
    !! \beta_{ijk} X_{i} X_{j} X_{k} + \sum_{i=1}^{n} \sum_{j=1 \\ j \neq i}^{n} 
    !! \beta_{iij} X_{i}^{2} X_{j} + \sum_{i=1}^{n} \beta_{iii} X_{i}^{3} $$
    !!
    !! The first line of terms (through the squared terms) is included for a 
    !! two-way model, and the remaining terms are added for a three-way model.  
    !! Each unique combination of factors appears exactly once.  For example, 
    !! with two factors, a two-way model is:
    !!
    !! $$ Y = \beta_{0} + \beta_{1} X_{1} + \beta_{2} X_{2} + \beta_{12} 
    !! X_{1} X_{2} + \beta_{11} X_{1}^{2} + \beta_{22} X_{2}^{2} $$
    !!
    !! The coefficients are ordered as: constant, main effects, cross terms 
    !! (i < j, lexicographic), squared terms, and (three-way only) triple 
    !! cross terms (i < j < k), X_i**2 * X_j terms (i, then j, with i /= j), 
    !! and cubed terms.
    integer, intent(in) :: nway
        !! The number of interaction levels.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! used to produce the results.
    real(real64), intent(in), dimension(:) :: y
        !! An M-element array containing the results from the M experiments.
    logical, intent(in), optional, target, dimension(:) :: map
        !! An optional array of the same size as beta that can be used to
        !! eliminate a parameter from the model (false), or keep a parameter
        !! in the model (true).  If not supplied, all parameters will be assumed
        !! to be part of the model as if the array were filled with all true
        !! values.
    real(real64), intent(in), optional :: alpha
        !! The significance level at which to evaluate the confidence 
        !! intervals.  The default value is 0.05 such that a 95% 
        !! confidence interval is calculated.
    type(doe_model) :: rst
        !! The resulting model.

    ! Local Variables
    integer :: i, j, m, n, nparam, nfactors
    logical, allocatable, target, dimension(:) :: nmap
    logical, pointer, dimension(:) :: mapptr
    real(real64) :: alph, nan
    real(real64), allocatable, dimension(:) :: coeffs, ymod, resid
    real(real64), allocatable, dimension(:,:) :: xc, c, cxt
    type(regression_statistics), allocatable, dimension(:) :: stats
    
    ! Initialization
    if (present(alpha)) then
        alph = alpha
    else
        alph = 5.0d-2
    end if
    m = size(x, 1)
    nfactors = size(x, 2)
    nan = ieee_value(nan, IEEE_QUIET_NAN)

    ! Input Checking
    if (nway < 1 .or. nway > 3) error stop FS_INVALID_INPUT_ERROR
    if (size(y) /= m) error stop FS_ARRAY_SIZE_ERROR

    ! Determine the parameter count
    nparam = doe_param_count(nway, nfactors)

    
    ! Set up the map parameters
    if (present(map)) then
        if (size(map) /= nparam) error stop FS_ARRAY_SIZE_ERROR
        mapptr => map
    else
        allocate(nmap(nparam), source = .true.)
        mapptr => nmap
    end if

    ! Update the parameter count
    n = nparam
    do i = 1, nparam
        if (.not.mapptr(i)) n = n - 1
    end do
    if (n < 1) then
        error stop FS_INVALID_INPUT_ERROR
    end if

    ! Local memory allocations
    allocate(xc(m, n), c(n, n), cxt(n, m), coeffs(n))

    ! Create the design matrix
    call doe_design_matrix(nway, x, mapptr, xc)

    ! Compute the covariance matrix
    c = covariance_matrix(xc)

    ! Solve the least-squares problem (N-by-1 result)
    call DGEMM("N", "T", n, m, n, 1.0d0, c, n, xc, m, 0.0d0, cxt, n) ! C * X**T
    call DGEMM("N", "N", n, 1, m, 1.0d0, cxt, n, y, m, 0.0d0, coeffs, n) ! (C * X**T) * Y

    ! Evaluate the model and compute the residuals
    ymod = matmul(xc, coeffs)
    resid = ymod - y

    ! Estimate parameter statistics
    stats = calculate_regression_statistics(resid, coeffs, c, alph)

    ! Update output
    rst%nway = nway
    rst%dof = m - n
    rst%mse = dot_product(resid, resid) / real(rst%dof, real64)
    allocate(rst%covariance, source = c)
    allocate(rst%coefficients(nparam))
    allocate(rst%stats(nparam))
    allocate(rst%map(nparam), source = mapptr)
    j = 0
    do i = 1, nparam
        if (mapptr(i)) then
            j = j + 1
            rst%coefficients(i) = coeffs(j)
            rst%stats(i) = stats(j)
        else
            rst%coefficients(i) = nan
            rst%stats(i)%confidence_interval = nan
            rst%stats(i)%probability = nan
            rst%stats(i)%standard_error = nan
            rst%stats(i)%t_statistic = nan
        end if
    end do
end function

! ------------------------------------------------------------------------------
function doe_model_covariance(mdl, x) result(rst)
    !! Computes the unscaled parameter covariance matrix 
    !! \( \left( X^{T} X \right)^{-1} \) for a model, where \( X \) is the 
    !! design matrix constructed from the supplied factor values.  Only the 
    !! parameters included in the model (map = true) are represented.
    class(doe_model), intent(in) :: mdl
        !! The DOE model defining the form of the design matrix.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! used to fit the model.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The P-by-P covariance matrix, where P is the number of parameters 
        !! included in the model.

    ! Local Variables
    integer :: m, p
    real(real64), allocatable, dimension(:,:) :: xc

    m = size(x, 1)
    p = count(mdl%map)
    if (p < 1) error stop FS_INVALID_INPUT_ERROR
    if (size(mdl%map) /= doe_param_count(mdl%nway, size(x, 2))) &
        error stop FS_ARRAY_SIZE_ERROR
    allocate(xc(m, p))
    call doe_design_matrix(mdl%nway, x, mdl%map, xc)
    rst = covariance_matrix(xc)
end function

! ------------------------------------------------------------------------------
subroutine doe_design_terms(mdl, x, terms)
    !! Evaluates the model terms (columns of the design matrix) included in
    !! a model at each of the supplied factor values.
    class(doe_model), intent(in) :: mdl
        !! The DOE model defining the terms.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors.
    real(real64), intent(out), dimension(:,:) :: terms
        !! The M-by-P matrix, where P is the number of parameters included in 
        !! the model, that will contain the terms.
    
    if (size(terms, 1) /= size(x, 1) .or. size(terms, 2) /= count(mdl%map)) &
        error stop FS_ARRAY_SIZE_ERROR
    call doe_design_matrix(mdl%nway, x, mdl%map, terms)
end subroutine

! ------------------------------------------------------------------------------
pure function doe_param_count(nway, n) result(rst)
    !! Computes the number of model parameters for an N-factor model.
    integer, intent(in) :: nway
        !! The number of interaction levels.
    integer, intent(in) :: n
        !! The number of factors.
    integer :: rst
        !! The number of parameters.

    rst = 1
    if (nway >= 1) rst = rst + n
    if (nway >= 2) rst = rst + n * (n + 1) / 2
    if (nway >= 3) rst = rst + n * (n + 1) * (n + 2) / 6
end function
! ------------------------------------------------------------------------------
subroutine doe_design_matrix(nway, x, map, c)
    !! This is an internal routine used to construct the design matrix for
    !! the DOE model of the following form:
    !!
    !! $$ Y = \beta_{0} + \sum_{i=1}^{n} \beta_{i} X_{i} + \sum_{i=1}^{n-1} 
    !! \sum_{j=i+1}^{n} \beta_{ij} X_{i} X_{j} + \sum_{i=1}^{n} \beta_{ii} 
    !! X_{i}^{2} + \sum_{i=1}^{n-2} \sum_{j=i+1}^{n-1} \sum_{k=j+1}^{n} 
    !! \beta_{ijk} X_{i} X_{j} X_{k} + \sum_{i=1}^{n} \sum_{j=1 \\ j \neq i}^{n} 
    !! \beta_{iij} X_{i}^{2} X_{j} + \sum_{i=1}^{n} \beta_{iii} X_{i}^{3} $$
    !!
    !! The first line of terms (through the squared terms) is included for a 
    !! two-way model, and the remaining terms are added for a three-way model.  
    !! Each unique combination of factors appears exactly once.  For example, 
    !! with two factors, a two-way model is:
    !!
    !! $$ Y = \beta_{0} + \beta_{1} X_{1} + \beta_{2} X_{2} + \beta_{12} 
    !! X_{1} X_{2} + \beta_{11} X_{1}^{2} + \beta_{22} X_{2}^{2} $$
    !!
    !! The coefficients are ordered as: constant, main effects, cross terms 
    !! (i < j, lexicographic), squared terms, and (three-way only) triple 
    !! cross terms (i < j < k), X_i**2 * X_j terms (i, then j, with i /= j), 
    !! and cubed terms.
    !!
    !! Up to a 3-way model is allowed.
    !!
    !! No error checking is provided.  It is assumed the arrays are sized 
    !! correctly.
    integer, intent(in) :: nway
    real(real64), intent(in), dimension(:,:) :: x
    logical, intent(in), dimension(:) :: map
    real(real64), intent(out), dimension(:,:) :: c

    ! Local Variables
    integer :: i, j, k, jj, kk, m, n, np

    ! Determine the number of model parameters
    np = 0
    do i = 1, size(map)
        if (map(i)) np = np + 1
    end do

    ! Additional Initialization
    m = size(x, 1)
    n = size(x, 2)

    ! DC Term
    if (map(1)) then
        c(:,1) = 1.0d0
        jj = 2
    else
        jj = 1
    end if

    ! Main Effect
    kk = 1
    if (nway >= 1) then
        do i = 1, n
            kk = kk + 1
            if (.not.map(kk)) cycle
            c(:,jj) = x(:,i)
            jj = jj + 1
        end do
    end if

    ! Two-Way (cross terms followed by squared terms)
    if (nway >= 2) then
        do i = 1, n - 1
            do j = i + 1, n
                kk = kk + 1
                if (.not.map(kk)) cycle
                c(:,jj) = x(:,i) * x(:,j)
                jj = jj + 1
            end do
        end do
        do i = 1, n
            kk = kk + 1
            if (.not.map(kk)) cycle
            c(:,jj) = x(:,i)**2
            jj = jj + 1
        end do
    end if

    ! Three-Way (triple cross terms, X_i**2 * X_j terms, and cubed terms)
    if (nway >= 3) then
        do i = 1, n - 2
            do j = i + 1, n - 1
                do k = j + 1, n
                    kk = kk + 1
                    if (.not.map(kk)) cycle
                    c(:,jj) = x(:,i) * x(:,j) * x(:,k)
                    jj = jj + 1
                end do
            end do
        end do
        do i = 1, n
            do j = 1, n
                if (i == j) cycle
                kk = kk + 1
                if (.not.map(kk)) cycle
                c(:,jj) = x(:,i)**2 * x(:,j)
                jj = jj + 1
            end do
        end do
        do i = 1, n
            kk = kk + 1
            if (.not.map(kk)) cycle
            c(:,jj) = x(:,i)**3
            jj = jj + 1
        end do
    end if
end subroutine

! ------------------------------------------------------------------------------
function doe_evaluate_model_1(nway, beta, x, map) result(rst)
    !! Evaluates the model of the following form.
    !!
    !! $$ Y = \beta_{0} + \sum_{i=1}^{n} \beta_{i} X_{i} + \sum_{i=1}^{n-1} 
    !! \sum_{j=i+1}^{n} \beta_{ij} X_{i} X_{j} + \sum_{i=1}^{n} \beta_{ii} 
    !! X_{i}^{2} + \sum_{i=1}^{n-2} \sum_{j=i+1}^{n-1} \sum_{k=j+1}^{n} 
    !! \beta_{ijk} X_{i} X_{j} X_{k} + \sum_{i=1}^{n} \sum_{j=1 \\ j \neq i}^{n} 
    !! \beta_{iij} X_{i}^{2} X_{j} + \sum_{i=1}^{n} \beta_{iii} X_{i}^{3} $$
    !!
    !! The first line of terms (through the squared terms) is included for a 
    !! two-way model, and the remaining terms are added for a three-way model.  
    !! Each unique combination of factors appears exactly once.  For example, 
    !! with two factors, a two-way model is:
    !!
    !! $$ Y = \beta_{0} + \beta_{1} X_{1} + \beta_{2} X_{2} + \beta_{12} 
    !! X_{1} X_{2} + \beta_{11} X_{1}^{2} + \beta_{22} X_{2}^{2} $$
    !!
    !! The coefficients are ordered as: constant, main effects, cross terms 
    !! (i < j, lexicographic), squared terms, and (three-way only) triple 
    !! cross terms (i < j < k), X_i**2 * X_j terms (i, then j, with i /= j), 
    !! and cubed terms.
    integer, intent(in) :: nway
        !! The number of interaction levels.  Currently, this algorithm supports
        !! a maximum of three-way interaction.
    real(real64), intent(in), dimension(:) :: beta
        !! The model coefficients.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! at which to evaluate the model.
    logical, intent(in), optional, target, dimension(:) :: map
        !! An optional array of the same size as beta that can be used to
        !! eliminate a parameter from the model (false), or keep a parameter
        !! in the model (true).  If not supplied, all parameters will be assumed
        !! to be part of the model as if the array were filled with all true
        !! values.
    real(real64), allocatable, dimension(:) :: rst
        !! The resulting M-element array.

    ! Local Variables
    integer :: m, n, nparam
    logical, pointer, dimension(:) :: mapptr
    logical, allocatable, target, dimension(:) :: nmap
    
    ! Initialization
    m = size(x, 1)
    n = size(x, 2)

    ! Input Checking
    if (nway < 1 .or. nway > 3) error stop FS_INVALID_INPUT_ERROR

    nparam = doe_param_count(nway, n)

    if (size(beta) /= nparam) error stop FS_ARRAY_SIZE_ERROR

    ! Memory Allocations
    allocate(rst(m))

    ! Set up the map parameters
    if (present(map)) then
        if (size(map) /= nparam) error stop FS_ARRAY_SIZE_ERROR
        mapptr => map
    else
        allocate(nmap(nparam), source = .true.)
        mapptr => nmap
    end if

    ! Process
    call doe_eval_engine(nway, beta, x, mapptr, rst)
end function

! ------------------------------------------------------------------------------
function doe_evaluate_model_2(mdl, x) result(rst)
    !! Evaluates the model of the following form.
    !!
    !! $$ Y = \beta_{0} + \sum_{i=1}^{n} \beta_{i} X_{i} + \sum_{i=1}^{n-1} 
    !! \sum_{j=i+1}^{n} \beta_{ij} X_{i} X_{j} + \sum_{i=1}^{n} \beta_{ii} 
    !! X_{i}^{2} + \sum_{i=1}^{n-2} \sum_{j=i+1}^{n-1} \sum_{k=j+1}^{n} 
    !! \beta_{ijk} X_{i} X_{j} X_{k} + \sum_{i=1}^{n} \sum_{j=1 \\ j \neq i}^{n} 
    !! \beta_{iij} X_{i}^{2} X_{j} + \sum_{i=1}^{n} \beta_{iii} X_{i}^{3} $$
    !!
    !! The first line of terms (through the squared terms) is included for a 
    !! two-way model, and the remaining terms are added for a three-way model.  
    !! Each unique combination of factors appears exactly once.  For example, 
    !! with two factors, a two-way model is:
    !!
    !! $$ Y = \beta_{0} + \beta_{1} X_{1} + \beta_{2} X_{2} + \beta_{12} 
    !! X_{1} X_{2} + \beta_{11} X_{1}^{2} + \beta_{22} X_{2}^{2} $$
    !!
    !! The coefficients are ordered as: constant, main effects, cross terms 
    !! (i < j, lexicographic), squared terms, and (three-way only) triple 
    !! cross terms (i < j < k), X_i**2 * X_j terms (i, then j, with i /= j), 
    !! and cubed terms.
    class(doe_model), intent(in) :: mdl
        !! The model to evaluate.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! at which to evaluate the model.
    real(real64), allocatable, dimension(:) :: rst
        !! The resulting M-element array.

    ! Process
    rst = doe_evaluate_model_1(mdl%nway, mdl%coefficients, x, mdl%map)
end function

! ------------------------------------------------------------------------------
subroutine doe_eval_engine(nway, beta, x, map, y)
    ! Driver routine for "doe_evaluate_model" that performs the actual 
    ! calculations but forgoes any error checking.  This should not be exposed 
    ! as part of the public API.
    integer, intent(in) :: nway
    real(real64), intent(in), dimension(:) :: beta
    real(real64), intent(in), dimension(:,:) :: x
    logical, intent(in), dimension(:) :: map
    real(real64), intent(out), dimension(:) :: y

    ! Local Variables
    integer :: i1, i2, n

    ! Initialization
    n = size(x, 2)
    if (map(1)) then
        y = beta(1)
    else
        y = 0.0d0
    end if

    ! Process
    if (nway >= 1) then
        i1 = 2
        i2 = i1 + n - 1
        call doe_eval_1(beta(i1:i2), x, map(i1:i2), y)
    end if
    if (nway >= 2) then
        i1 = i2 + 1
        i2 = i1 + n * (n + 1) / 2 - 1
        call doe_eval_2(beta(i1:i2), x, map(i1:i2), y)
    end if
    if (nway >= 3) then
        i1 = i2 + 1
        i2 = i1 + n * (n + 1) * (n + 2) / 6 - 1
        call doe_eval_3(beta(i1:i2), x, map(i1:i2), y)
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine doe_eval_1(beta, x, map, y)
    !! Evaluates the main effect term.
    !!
    !! $$ Y = Y + \sum_{i=1}^{n} \beta_{i} X_{i} $$
    real(real64), intent(in), dimension(:) :: beta
        !! The model coefficients for just this portion of the model.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! at which to evaluate the model.
    logical, intent(in), dimension(:) :: map
        !! The usage map corresponding to the model coefficients for just this
        !! portion of the model.
    real(real64), intent(inout), dimension(:) :: y
        !! On input, an M-element array containing the existing portion of the 
        !! model.  On output, this array is updated to include the main effects.

    ! Local Variables
    integer :: i, n

    ! Initialization
    n = size(x, 2)

    ! Process
    do i = 1, n
        if (.not.map(i)) cycle
        y = y + beta(i) * x(:,i)
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine doe_eval_2(beta, x, map, y)
    !! Evaluates the two-way interaction term.
    !!
    !! $ Y = Y + \sum_{i=1}^{n-1} \sum_{j=i+1}^{n} \beta_{ij} X_{i} X_{j} +`r`n    !! \sum_{i=1}^{n} \beta_{ii} X_{i}^{2} $
    real(real64), intent(in), dimension(:) :: beta
        !! The model coefficients for just this portion of the model.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! at which to evaluate the model.
    logical, intent(in), dimension(:) :: map
        !! The usage map corresponding to the model coefficients for just this
        !! portion of the model.
    real(real64), intent(inout), dimension(:) :: y
        !! On input, an M-element array containing the existing portion of the 
        !! model.  On output, this array is updated to include the two-way
        !! interactions.

    ! Local Variables
    integer :: i, j, k, n

    ! Initialization
    n = size(x, 2)

    ! Process
    k = 0
    do i = 1, n - 1
        do j = i + 1, n
            k = k + 1
            if (.not.map(k)) cycle
            y = y + beta(k) * x(:,i) * x(:,j)
        end do
    end do
    do i = 1, n
        k = k + 1
        if (.not.map(k)) cycle
        y = y + beta(k) * x(:,i)**2
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine doe_eval_3(beta, x, map, y)
    !! Evaluates the three-way interaction term.
    !!
    !! $ Y = Y + \sum_{i=1}^{n-2} \sum_{j=i+1}^{n-1} \sum_{k=j+1}^{n}`r`n    !! \beta_{ijk} X_{i} X_{j} X_{k} + \sum_{i=1}^{n} \sum_{j=1 \\ j \neq i}^{n}`r`n    !! \beta_{iij} X_{i}^{2} X_{j} + \sum_{i=1}^{n} \beta_{iii} X_{i}^{3} $
    real(real64), intent(in), dimension(:) :: beta
        !! The model coefficients for just this portion of the model.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! at which to evaluate the model.
    logical, intent(in), dimension(:) :: map
        !! The usage map corresponding to the model coefficients for just this
        !! portion of the model.
    real(real64), intent(inout), dimension(:) :: y
        !! On input, an M-element array containing the existing portion of the 
        !! model.  On output, this array is updated to include the three-way
        !! interactions.

    ! Local Variables
    integer :: i, j, k, ii, n

    ! Initialization
    n = size(x, 2)

    ! Process
    ii = 0
    do i = 1, n - 2
        do j = i + 1, n - 1
            do k = j + 1, n
                ii = ii + 1
                if (.not.map(ii)) cycle
                y = y + beta(ii) * x(:,i) * x(:,j) * x(:,k)
            end do
        end do
    end do
    do i = 1, n
        do j = 1, n
            if (i == j) cycle
            ii = ii + 1
            if (.not.map(ii)) cycle
            y = y + beta(ii) * x(:,i)**2 * x(:,j)
        end do
    end do
    do i = 1, n
        ii = ii + 1
        if (.not.map(ii)) cycle
        y = y + beta(ii) * x(:,i)**3
    end do
end subroutine

end module fstats_doe_models

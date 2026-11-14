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

module fstats_doe_prediction
    !! DOE predictions with confidence and prediction intervals.
    !!
    !! This module evaluates a doe_model at a set of factor values and returns
    !! the predicted responses along with confidence and prediction intervals.
    !! With \(\mathbf{x}_{0}\) the design-matrix row (model terms) for a 
    !! prediction point, \(C = \left( X^{T} X \right)^{-1}\) the unscaled 
    !! parameter covariance of the fit, \(MSE\) the residual mean squared 
    !! error, and \(t\) the critical value of Student's t distribution with
    !! \(\nu = m - p\) degrees of freedom at probability \(1 - \alpha / 2\), 
    !! the intervals are
    !!
    !! $$ \hat{y} \pm t \sqrt{MSE \cdot \mathbf{x}_{0}^{T} C \mathbf{x}_{0}} $$
    !! $$ \hat{y} \pm t \sqrt{MSE \left( 1 + \mathbf{x}_{0}^{T} C 
    !! \mathbf{x}_{0} \right)} $$
    !!
    !! for the mean response (confidence interval) and for a new observation
    !! (prediction interval), respectively.  The factor values must be in the 
    !! same scale (e.g. coded on [-1, 1]) used when the model was fit.
    !!
    !! The covariance matrix and the residual statistics are stored in the 
    !! model by doe_fit_model.  For a model that was not created by 
    !! doe_fit_model, supply the factor values used to fit the model (x_fit) 
    !! and the residual mean squared error (residual_mse) to the routines.
    use iso_fortran_env, only : int32, real64
    use fstats_errors, only : FS_INVALID_INPUT_ERROR
    use fstats_experimental_design, only : doe_model, doe_prediction
    use fstats_doe_models, only : doe_evaluate_model, doe_model_covariance, &
        doe_design_terms
    use fstats_t_distribution, only : t_distribution
    implicit none
    private

    public :: doe_predict
    public :: doe_predict_enhanced
contains
! ------------------------------------------------------------------------------
function doe_predict(mdl, x, alpha, x_fit, residual_mse) result(pred)
    !! Computes predictions with confidence and prediction intervals.
    !!
    !! The intervals are computed as described in the module documentation 
    !! using the leverage of each prediction point and the critical value of
    !! Student's t distribution.  The covariance matrix and the residual mean
    !! squared error are taken from the model if it was created by 
    !! doe_fit_model; otherwise, x_fit and residual_mse must be supplied.  If
    !! either is supplied, it takes precedence over the value stored in the 
    !! model.  The routine halts with an error if the required information is
    !! not available, or if the residual degrees of freedom are not positive.
    class(doe_model), intent(in) :: mdl
        !! The DOE model.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! at which to evaluate the model.
    real(real64), intent(in), optional :: alpha
        !! The significance level, in the range (0, 1).  The default is 0.05
        !! (95% intervals).
    real(real64), intent(in), optional, dimension(:,:) :: x_fit
        !! An optional K-by-N matrix containing the factor values used to fit
        !! the model.  This is used to compute the parameter covariance and 
        !! the residual degrees of freedom \(K - p\).
    real(real64), intent(in), optional :: residual_mse
        !! An optional residual mean squared error of the fit.
    type(doe_prediction) :: pred
        !! A doe_prediction object containing the M predicted values, the 
        !! lower and upper confidence bounds, the lower and upper prediction 
        !! bounds, and the confidence level.

    ! Local Variables
    integer :: i, m, p, dof
    real(real64) :: alph, mse, tcrit, q
    real(real64), allocatable, dimension(:,:) :: c, xt, w
    type(t_distribution) :: dist

    ! Initialization
    m = size(x, 1)
    p = count(mdl%map)
    alph = 5.0d-2
    if (present(alpha)) alph = alpha
    if (alph <= 0.0d0 .or. alph >= 1.0d0) error stop FS_INVALID_INPUT_ERROR

    ! Obtain the covariance matrix and residual statistics
    if (present(x_fit)) then
        c = doe_model_covariance(mdl, x_fit)
        dof = size(x_fit, 1) - p
    else if (allocated(mdl%covariance) .and. mdl%dof > 0) then
        c = mdl%covariance
        dof = mdl%dof
    else
        error stop FS_INVALID_INPUT_ERROR
    end if
    if (dof <= 0) error stop FS_INVALID_INPUT_ERROR
    if (present(residual_mse)) then
        mse = residual_mse
    else if (mdl%dof > 0) then
        mse = mdl%mse
    else
        error stop FS_INVALID_INPUT_ERROR
    end if
    if (mse < 0.0d0) error stop FS_INVALID_INPUT_ERROR

    ! Critical value
    dist%dof = real(dof, real64)
    tcrit = dist%standardized_variable(1.0d0 - 0.5d0 * alph)

    ! Predictions and leverage
    allocate(xt(m, p), w(m, p))
    call doe_design_terms(mdl, x, xt)
    w = matmul(xt, c)
    allocate(pred%predicted_values, source = doe_evaluate_model(mdl, x))
    allocate(pred%confidence_lower(m), pred%confidence_upper(m))
    allocate(pred%prediction_lower(m), pred%prediction_upper(m))
    do i = 1, m
        q = dot_product(w(i,:), xt(i,:))
        pred%confidence_lower(i) = pred%predicted_values(i) &
            - tcrit * sqrt(mse * q)
        pred%confidence_upper(i) = pred%predicted_values(i) &
            + tcrit * sqrt(mse * q)
        pred%prediction_lower(i) = pred%predicted_values(i) &
            - tcrit * sqrt(mse * (1.0d0 + q))
        pred%prediction_upper(i) = pred%predicted_values(i) &
            + tcrit * sqrt(mse * (1.0d0 + q))
    end do
    pred%confidence_level = 1.0d0 - alph
end function

! ------------------------------------------------------------------------------
function doe_predict_enhanced(mdl, x, alpha, residual_mse, x_fit) result(pred)
    !! Computes predictions with confidence and prediction intervals using a 
    !! user-supplied residual mean squared error.  This is equivalent to
    !! doe_predict, and is retained for convenience when the residual mean 
    !! squared error is computed by the caller.
    class(doe_model), intent(in) :: mdl
        !! The DOE model.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! at which to evaluate the model.
    real(real64), intent(in), optional :: alpha
        !! The significance level, in the range (0, 1).  The default is 0.05
        !! (95% intervals).
    real(real64), intent(in), optional :: residual_mse
        !! An optional residual mean squared error of the fit.  If omitted, 
        !! the value stored in the model by doe_fit_model is used.
    real(real64), intent(in), optional, dimension(:,:) :: x_fit
        !! An optional K-by-N matrix containing the factor values used to fit
        !! the model.  This is required if the model was not created by 
        !! doe_fit_model.
    type(doe_prediction) :: pred
        !! A doe_prediction object containing the M predicted values, the 
        !! lower and upper confidence bounds, the lower and upper prediction 
        !! bounds, and the confidence level.

    pred = doe_predict(mdl, x, alpha, x_fit, residual_mse)
end function

end module fstats_doe_prediction

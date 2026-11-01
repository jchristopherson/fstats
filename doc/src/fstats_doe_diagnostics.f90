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

module fstats_doe_diagnostics
    !! DOE goodness of fit, residual analysis, comparison, and ANOVA.
    !!
    !! The routines in this module operate on a fitted (or user-defined) 
    !! doe_model and the data used to evaluate it.  In each case, the model is
    !! evaluated at the M-by-N factor matrix X using doe_evaluate_model, and 
    !! the result is compared with the M observed responses Y.  The factor 
    !! values must be supplied in the same scale (e.g. coded on [-1, 1]) used 
    !! when the model was fit.
    !!
    !! The number of model parameters, \(p\), is taken as the number of 
    !! coefficients enabled in the model's usage map (i.e. `count(mdl%map)`).
    use iso_fortran_env, only : int32, real64
    use fstats_experimental_design, only : doe_model, doe_diagnostics, &
        doe_residuals, doe_comparison_result, doe_anova_table
    use fstats_doe_models, only : doe_evaluate_model
    implicit none
    private

    public :: doe_model_diagnostics
    public :: doe_residuals_analysis
    public :: doe_compare_models
    public :: doe_model_anova
contains
! ------------------------------------------------------------------------------
function doe_model_diagnostics(mdl, x, y) result(diag)
    !! Computes model diagnostics and goodness-of-fit metrics.
    !!
    !! With \(m\) observations, \(p\) model parameters, observed responses 
    !! \(y_i\), model predictions \(\hat{y}_i\), and mean response 
    !! \(\bar{y}\), the following quantities are computed.
    !!
    !! $$ SS_{res} = \sum_{i=1}^{m} (y_i - \hat{y}_i)^{2} $$
    !! $$ SS_{tot} = \sum_{i=1}^{m} (y_i - \bar{y})^{2} $$
    !! $$ R^{2} = 1 - \frac{SS_{res}}{SS_{tot}} $$
    !! $$ R^{2}_{adj} = 1 - \frac{SS_{res} / (m - p)}{SS_{tot} / (m - 1)} $$
    !! $$ RMSE = \sqrt{\frac{SS_{res}}{m}} $$
    !! $$ s = \sqrt{\frac{SS_{res}}{m - p}} $$
    !! $$ F = \frac{(SS_{tot} - SS_{res}) / (p - 1)}{SS_{res} / (m - p)} $$
    !!
    !! Degenerate cases are handled as follows: \(R^{2}\) is zero if 
    !! \(SS_{tot} = 0\); the adjusted \(R^{2}\) and residual standard error 
    !! are zero if \(m \le p\); and the F-statistic is zero (with a p-value of
    !! 1) if \(p = 1\) or \(m \le p\).
    !!
    !! Note, the F-statistic p-value is not currently computed.  When the 
    !! F-statistic is calculated, the p-value is reported as zero.
    class(doe_model), intent(in) :: mdl
        !! The fitted DOE model.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! used in model fitting.
    real(real64), intent(in), dimension(:) :: y
        !! The M-element array containing the observed responses.  This array
        !! must have the same number of rows as x.
    type(doe_diagnostics) :: diag
        !! A doe_diagnostics object containing the R-squared, adjusted
        !! R-squared, RMSE, residual standard error, F-statistic, mean
        !! response, observation count, and parameter count.

    ! Local Variables
    integer(int32) :: m, p, i
    real(real64) :: ss_total, ss_residual, y_pred, y_mean, mse
    real(real64), allocatable :: residuals(:)
    real(real64) :: f_stat, ss_model

    ! Initialization
    m = size(y)
    p = count(mdl%map)
    if (p < 1) p = 1  ! At least intercept

    y_mean = sum(y) / real(m, real64)
    
    ! Calculate residuals
    allocate(residuals(m))
    residuals = y - doe_evaluate_model(mdl, x)

    ! Calculate sum of squares
    ss_residual = sum(residuals**2)
    ss_total = sum((y - y_mean)**2)
    ss_model = ss_total - ss_residual

    ! R-squared
    if (ss_total > 0.0d0) then
        diag%r_squared = 1.0d0 - (ss_residual / ss_total)
    else
        diag%r_squared = 0.0d0
    end if

    ! Adjusted R-squared
    if (m - p > 0) then
        diag%r_squared_adjusted = 1.0d0 - &
            (ss_residual / real(m - p, real64)) / &
            (ss_total / real(m - 1, real64))
    else
        diag%r_squared_adjusted = 0.0d0
    end if

    ! RMSE and residual standard error
    if (m > 0) then
        diag%rmse = sqrt(ss_residual / real(m, real64))
    end if
    if (m - p > 0) then
        mse = ss_residual / real(m - p, real64)
        diag%residual_std_error = sqrt(mse)
    else
        diag%residual_std_error = 0.0d0
    end if

    ! F-statistic and p-value
    if (p > 1 .and. m - p > 0) then
        f_stat = (ss_model / real(p - 1, real64)) / &
                 (ss_residual / real(m - p, real64))
        diag%f_statistic = f_stat
        ! P-value calculation would require F-distribution CDF (not implemented here)
        diag%f_p_value = 0.0d0
    else
        diag%f_statistic = 0.0d0
        diag%f_p_value = 1.0d0
    end if

    ! Additional info
    diag%mean_response = y_mean
    diag%n_observations = m
    diag%n_parameters = p
end function

! ------------------------------------------------------------------------------
function doe_residuals_analysis(mdl, x, y) result(resid)
    !! Computes residual analysis data for model diagnostics.
    !!
    !! The residuals are defined as \(e_i = y_i - \hat{y}_i\).  The mean and 
    !! sample standard deviation (using \(m - 1\) in the denominator) of the 
    !! residuals are computed, and the standardized residuals are then 
    !! \(z_i = (e_i - \bar{e}) / s_e\).  If the standard deviation is zero 
    !! (or only a single observation is provided), the standard deviation and
    !! standardized residuals are reported as zero.
    class(doe_model), intent(in) :: mdl
        !! The fitted DOE model.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! at which to evaluate the model.
    real(real64), intent(in), dimension(:) :: y
        !! The M-element array containing the observed responses.  This array
        !! must have the same number of rows as x.
    type(doe_residuals) :: resid
        !! A doe_residuals object containing the observed values, predicted
        !! values, residuals, standardized residuals, and the residual mean
        !! and standard deviation.

    ! Local Variables
    integer(int32) :: m, i
    real(real64), allocatable :: y_pred(:)

    ! Initialization
    m = size(y)
    allocate(resid%observed_values, source=y)

    ! Calculate predictions
    y_pred = doe_evaluate_model(mdl, x)
    allocate(resid%predicted_values, source=y_pred)

    ! Calculate residuals
    allocate(resid%residuals(m))
    resid%residuals = y - y_pred

    ! Mean and std of residuals
    resid%residual_mean = sum(resid%residuals) / real(m, real64)
    if (m > 1) then
        resid%residual_std = sqrt(sum((resid%residuals - resid%residual_mean)**2) / &
                                  real(m - 1, real64))
    else
        resid%residual_std = 0.0d0
    end if

    ! Standardized residuals
    allocate(resid%standardized_residuals(m))
    if (resid%residual_std > 0.0d0) then
        resid%standardized_residuals = (resid%residuals - resid%residual_mean) / &
                                       resid%residual_std
    else
        resid%standardized_residuals = 0.0d0
    end if

end function

! ------------------------------------------------------------------------------
function doe_compare_models(mdl1, mdl2, x, y) result(comp)
    !! Compares two DOE models using the same data set.
    !!
    !! The residual sum of squares (RSS) of each model is computed by 
    !! evaluating each model at x and comparing the result to y.  Model 1 is
    !! treated as the "full" model, and model 2 as the "reduced" model.
    !!
    !! NOTE: This routine is currently a placeholder.  Only the two RSS values
    !! are computed from the data.  The remaining outputs are fixed values 
    !! that do not depend on the models or the data: the F-statistic is 1.5, 
    !! the p-value is 0.1, no significant difference is reported, the degrees 
    !! of freedom are M - 3 (full), M - 2 (reduced), and 1 (difference), and 
    !! the conclusion text is constant.  A proper nested model F-test is
    !!
    !! $$ F = \frac{(RSS_{reduced} - RSS_{full}) / (df_{reduced} - 
    !! df_{full})}{RSS_{full} / df_{full}} $$
    class(doe_model), intent(in) :: mdl1
        !! The full model.
    class(doe_model), intent(in) :: mdl2
        !! The reduced model.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! at which to evaluate both models.
    real(real64), intent(in), dimension(:) :: y
        !! The M-element array containing the observed responses.  This array
        !! must have the same number of rows as x.
    type(doe_comparison_result) :: comp
        !! A doe_comparison_result object containing the comparison results.

    comp%f_statistic = 1.5d0
    comp%p_value = 0.1d0
    comp%significant_difference = .false.
    comp%rss_full = sum((y - doe_evaluate_model(mdl1, x))**2)
    comp%rss_reduced = sum((y - doe_evaluate_model(mdl2, x))**2)
    comp%df_full = size(y) - 3
    comp%df_reduced = size(y) - 2
    comp%df_diff = 1
    comp%conclusion = "Models similar; use simpler one."
end function

! ------------------------------------------------------------------------------
function doe_model_anova(mdl, x, y) result(anova)
    !! Computes an analysis of variance (ANOVA) table for a DOE model.
    !!
    !! The sums of squares are computed from the data as
    !!
    !! $$ SS_{tot} = \sum_{i=1}^{m} (y_i - \bar{y})^{2}, \quad 
    !! SS_{res} = \sum_{i=1}^{m} (y_i - \hat{y}_i)^{2}, \quad 
    !! SS_{model} = SS_{tot} - SS_{res} $$
    !!
    !! with \(R^{2} = 1 - SS_{res} / SS_{tot}\), \(MS = SS / df\), and 
    !! \(F = MS_{model} / MS_{res}\).
    !!
    !! NOTE: This routine is currently only partially complete.  The degrees of
    !! freedom are fixed at 2 for the model and M - 3 for the residual (rather
    !! than being derived from the number of model parameters), and the 
    !! p-value is fixed at 0.01 rather than being computed from the 
    !! F-distribution.  The total degrees of freedom is M - 1.
    class(doe_model), intent(in) :: mdl
        !! The DOE model to analyze.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix containing the M values of each of the N factors
        !! at which to evaluate the model.
    real(real64), intent(in), dimension(:) :: y
        !! The M-element array containing the observed responses.  This array
        !! must have the same number of rows as x.
    type(doe_anova_table) :: anova
        !! A doe_anova_table object containing the ANOVA table.

    ! Local Variables
    real(real64), allocatable :: y_pred(:)

    y_pred = doe_evaluate_model(mdl, x)
    
    anova%ss_total = sum((y - sum(y)/real(size(y), real64))**2)
    anova%ss_residual = sum((y - y_pred)**2)
    anova%ss_model = anova%ss_total - anova%ss_residual
    anova%df_total = size(y) - 1
    anova%df_model = 2
    anova%df_residual = size(y) - 3
    anova%ms_model = anova%ss_model / real(max(1, anova%df_model), real64)
    anova%ms_residual = anova%ss_residual / real(max(1, anova%df_residual), real64)
    anova%f_statistic = anova%ms_model / max(anova%ms_residual, 1.0d-10)
    anova%p_value = 0.01d0
    anova%r_squared = 1.0d0 - (anova%ss_residual / max(anova%ss_total, 1.0d-10))
end function

end module fstats_doe_diagnostics

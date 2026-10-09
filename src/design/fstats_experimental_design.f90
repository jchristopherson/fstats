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

module fstats_experimental_design
    !! Shared types for experimental design generation, modeling, and analysis.
    use iso_fortran_env
    use fstats_regression, only : regression_statistics
    implicit none
    private
    
    public :: doe_model
    public :: doe_diagnostics
    public :: doe_prediction
    public :: doe_residuals
    public :: doe_efficiency_metrics
    public :: doe_rsm_model
    public :: doe_optimization_result
    public :: doe_comparison_result
    public :: doe_anova_table

    type doe_model
        !! A model used to represent a design of experiments result.  The model
        !! is of the following form.
        !!
        !! $$ Y = \beta_{0} + \sum_{i=1}^{n} \beta_{i} X_{i} + \sum_{i=1}^{n} 
        !! \sum_{j=1 \\ i \neq j}^{n} \beta_{ij} X_{i} X_{j} + \sum_{i=1}^{n} 
        !! \sum_{j=1}^{n} \sum_{k=1 \\ i \neq j \neq k}^{n} \beta_{ijk} X_{i} 
        !! X_{j} X_{k} + ... $$
        integer(int32) :: nway
            !! The number of interaction levels.
        real(real64), allocatable, dimension(:) :: coefficients
            !! The model coefficients.
        type(regression_statistics), allocatable, dimension(:) :: stats
            !! Statistical information for each model parameter.
        logical, allocatable, dimension(:) :: map
            !! An array denoting if a model coefficient should be included
            !! as part of the model (true), or neglected (false).
    end type

    type doe_diagnostics
        !! Model diagnostics and goodness-of-fit metrics.
        real(real64) :: r_squared
            !! The coefficient of determination (RÂ²), range [0, 1].
        real(real64) :: r_squared_adjusted
            !! The adjusted RÂ² accounting for model complexity.
        real(real64) :: rmse
            !! Root mean square error.
        real(real64) :: residual_std_error
            !! Residual standard error (standard deviation of residuals).
        real(real64) :: f_statistic
            !! Overall F-statistic for the model.
        real(real64) :: f_p_value
            !! P-value for the overall F-statistic.
        real(real64) :: mean_response
            !! Mean of the response variable.
        integer(int32) :: n_observations
            !! Number of observations.
        integer(int32) :: n_parameters
            !! Number of model parameters (including intercept).
    end type

    type doe_prediction
        !! Prediction with uncertainty quantification.
        real(real64), allocatable, dimension(:) :: predicted_values
            !! Predicted response values.
        real(real64), allocatable, dimension(:) :: confidence_lower
            !! Lower confidence interval bounds.
        real(real64), allocatable, dimension(:) :: confidence_upper
            !! Upper confidence interval bounds.
        real(real64), allocatable, dimension(:) :: prediction_lower
            !! Lower prediction interval bounds.
        real(real64), allocatable, dimension(:) :: prediction_upper
            !! Upper prediction interval bounds.
        real(real64) :: confidence_level
            !! Confidence level (e.g., 0.95 for 95% CI).
    end type

    type doe_residuals
        !! Residual analysis data.
        real(real64), allocatable, dimension(:) :: residuals
            !! Raw residuals: y - y_predicted.
        real(real64), allocatable, dimension(:) :: standardized_residuals
            !! Standardized residuals for outlier detection.
        real(real64), allocatable, dimension(:) :: predicted_values
            !! Predicted response values.
        real(real64), allocatable, dimension(:) :: observed_values
            !! Observed response values.
        real(real64) :: residual_mean
            !! Mean of residuals (should be ~0).
        real(real64) :: residual_std
            !! Standard deviation of residuals.
    end type

    type doe_efficiency_metrics
        !! Design efficiency metrics for evaluating design quality.
        real(real64) :: d_efficiency
            !! D-efficiency: \((|X'X|^(1/p))^(1/n)\) where p=params, n=runs.
            !! Range [0,1]. Higher is better (max=1 for orthogonal designs).
        real(real64) :: a_efficiency
            !! A-efficiency: \(p / trace((X'X)^-1)\). Range [0,1]. 
            !! Higher is better.
        real(real64) :: g_efficiency
            !! G-efficiency: 1 - (max_prediction_variance / avg_prediction_variance)
        real(real64) :: orthogonality
            !! Orthogonality measure: \(1.0\) if perfectly orthogonal, \(<1.0\) 
            !! otherwise.
        logical :: is_orthogonal
            !! True if design is perfectly orthogonal.
        integer(int32) :: n_runs
            !! Number of design runs.
        integer(int32) :: n_factors
            !! Number of factors.
        integer(int32) :: n_parameters
            !! Number of model parameters.
    end type

    type doe_rsm_model
        !! Response Surface Model (quadratic model for RSM).
        type(doe_model) :: base_model
            !! Base fitted model.
        real(real64), allocatable, dimension(:) :: linear_coeff
            !! Linear coefficients for each factor.
        real(real64), allocatable, dimension(:) :: quadratic_coeff
            !! Quadratic coefficients (main effects squared).
        real(real64), allocatable, dimension(:) :: interaction_coeff
            !! Interaction coefficients.
        real(real64) :: intercept
            !! Model intercept.
        real(real64) :: response_at_center
            !! Predicted response at design center (coded 0,0,...,0).
        integer(int32) :: n_factors
            !! Number of factors.
    end type

    type doe_optimization_result
        !! Results from RSM-based optimization.
        real(real64), allocatable, dimension(:) :: optimal_coded_factors
            !! Optimal factor settings (in coded scale).
        real(real64), allocatable, dimension(:) :: optimal_natural_factors
            !! Optimal factor settings (in natural scale).
        real(real64) :: optimal_response
            !! Predicted response at optimal point.
        integer(int32) :: iteration_count
            !! Number of iterations to converge.
        logical :: converged
            !! Whether optimization converged.
        character(len=256) :: method
            !! Optimization method used.
        real(real64) :: convergence_tolerance
            !! Tolerance used for convergence.
    end type

    type doe_comparison_result
        !! Results from comparing two models.
        real(real64) :: f_statistic
            !! F-statistic for model comparison.
        real(real64) :: p_value
            !! P-value for the F-test.
        real(real64) :: rss_full
            !! Residual sum of squares for full model.
        real(real64) :: rss_reduced
            !! Residual sum of squares for reduced model.
        integer(int32) :: df_full
            !! Degrees of freedom for full model.
        integer(int32) :: df_reduced
            !! Degrees of freedom for reduced model.
        integer(int32) :: df_diff
            !! Difference in degrees of freedom.
        logical :: significant_difference
            !! True if models differ significantly (p < 0.05).
        character(len=256) :: conclusion
            !! Interpretation of comparison results.
    end type

    type doe_anova_table
        !! ANOVA table for overall model fit assessment.
        real(real64) :: ss_total
            !! Total sum of squares.
        real(real64) :: ss_model
            !! Model sum of squares.
        real(real64) :: ss_residual
            !! Residual sum of squares.
        integer(int32) :: df_total
            !! Total degrees of freedom.
        integer(int32) :: df_model
            !! Model degrees of freedom.
        integer(int32) :: df_residual
            !! Residual degrees of freedom.
        real(real64) :: ms_model
            !! Model mean square.
        real(real64) :: ms_residual
            !! Residual mean square.
        real(real64) :: f_statistic
            !! F-statistic.
        real(real64) :: p_value
            !! P-value for the F-test.
        real(real64) :: r_squared
            !! R-squared value.
    end type
end module fstats_experimental_design

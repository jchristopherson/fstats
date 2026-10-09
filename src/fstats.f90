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

module fstats
    !! FSTATS is a modern Fortran statistical library containing routines for 
    !! computing basic statistical properties, hypothesis testing, regression, 
    !! special functions, and experimental design.
    use iso_fortran_env
    use fstats_special_functions
    use fstats_descriptive_statistics
    use fstats_hypothesis
    use fstats_distributions
    use fstats_t_distribution
    use fstats_normal_distribution
    use fstats_f_distribution
    use fstats_chi_squared_distribution
    use fstats_binomial_distribution
    use fstats_log_normal_distribution
    use fstats_poisson_distribution
    use fstats_multivariate_normal_distribution
    use fstats_anova
    use fstats_helper_routines
    use fstats_regression
    use fstats_linear_regression
    use fstats_levenberg_marquardt
    use fstats_experimental_design
    use fstats_doe_designs
    use fstats_doe_coding
    use fstats_doe_models
    use fstats_doe_diagnostics
    use fstats_doe_prediction
    use fstats_doe_response_surface
    use fstats_allan
    use fstats_bootstrap
    use fstats_sampling
    use fstats_smoothing
    use fstats_mcmc
    use fstats_interp
    use fstats_missing_data
    use fstats_msa
    use fstats_robust_statistics
end module
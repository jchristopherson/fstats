program doe_model_comparison_example
    !! Example demonstrating model comparison using ANOVA and F-tests
    use iso_fortran_env
    use fstats
    implicit none

    integer(int32) :: i, j, k
    integer(int32), parameter :: npts = 9
    real(real64) :: x(npts, 2), y(npts)
    real(real64) :: beta_simple(3), beta_full(5)
    type(doe_model) :: mdl_simple, mdl_full
    type(doe_comparison_result) :: comp
    type(doe_anova_table) :: tbl_simple, tbl_full

    print '(A)', "================================================"
    print '(A)', "Design of Experiments - Model Comparison"
    print '(A)', "================================================"
    print '(A)', ""

    ! Generate 3x3 factorial design
    k = 1
    do i = 1, 3
        do j = 1, 3
            x(k, 1) = real(i - 1, real64)
            x(k, 2) = real(j - 1, real64)
            y(k) = 1.0d0 + 2.0d0*x(k,1) + 1.5d0*x(k,2) + 0.3d0*x(k,1)*x(k,2)
            k = k + 1
        end do
    end do

    ! =====================================================================
    ! Model 1: Simple Model (Main Effects Only)
    ! =====================================================================
    beta_simple = [1.0d0, 2.0d0, 1.5d0]
    allocate(mdl_simple%coefficients(3), source=beta_simple)
    allocate(mdl_simple%stats(3))
    allocate(mdl_simple%map(3), source=[.true., .true., .true.])
    mdl_simple%nway = 1

    print '(A)', "MODEL 1: Simple Model (Main Effects Only)"
    print '(A)', "-" // repeat("-", 50)
    print '(A)', "  Equation: y = beta0 + beta1*x1 + beta2*x2"
    print '(A,F8.4,A,F8.4,A,F8.4)', "  Coefficients: b0=", beta_simple(1), &
        ", b1=", beta_simple(2), ", b2=", beta_simple(3)
    print '(A)', ""

    tbl_simple = doe_model_anova(mdl_simple, x, y)
    print '(A)', "  ANOVA Results:"
    print '(A,F12.4)', "    SS_total     = ", tbl_simple%ss_total
    print '(A,F12.4)', "    SS_model     = ", tbl_simple%ss_model
    print '(A,F12.4)', "    SS_residual  = ", tbl_simple%ss_residual
    print '(A,I0)', "    df_model     = ", tbl_simple%df_model
    print '(A,I0)', "    df_residual  = ", tbl_simple%df_residual
    print '(A,F12.6)', "    F-statistic  = ", tbl_simple%f_statistic
    print '(A,F8.6)', "    R-squared    = ", tbl_simple%r_squared
    print '(A)', ""

    ! =====================================================================
    ! Model 2: Full Model (With Interaction)
    ! =====================================================================
    beta_full = [1.0d0, 2.0d0, 1.5d0, 0.3d0, 0.0d0]
    allocate(mdl_full%coefficients(5), source=beta_full)
    allocate(mdl_full%stats(5))
    allocate(mdl_full%map(5), source=[.true., .true., .true., .true., .false.])
    mdl_full%nway = 2

    print '(A)', "MODEL 2: Full Model (With Interaction)"
    print '(A)', "-" // repeat("-", 50)
    print '(A)', "  Equation: y = beta0 + beta1*x1 + beta2*x2 + beta3*x1*x2"
    print '(A,F8.4,A,F8.4,A,F8.4,A,F8.4)', "  Coefficients: b0=", beta_full(1), &
        ", b1=", beta_full(2), ", b2=", beta_full(3), ", b3=", beta_full(4)
    print '(A)', ""

    tbl_full = doe_model_anova(mdl_full, x, y)
    print '(A)', "  ANOVA Results:"
    print '(A,F12.4)', "    SS_total     = ", tbl_full%ss_total
    print '(A,F12.4)', "    SS_model     = ", tbl_full%ss_model
    print '(A,F12.4)', "    SS_residual  = ", tbl_full%ss_residual
    print '(A,I0)', "    df_model     = ", tbl_full%df_model
    print '(A,I0)', "    df_residual  = ", tbl_full%df_residual
    print '(A,F12.6)', "    F-statistic  = ", tbl_full%f_statistic
    print '(A,F8.6)', "    R-squared    = ", tbl_full%r_squared
    print '(A)', ""

    ! =====================================================================
    ! Compare Models Using F-Test
    ! =====================================================================
    comp = doe_compare_models(mdl_simple, mdl_full, x, y)

    print '(A)', "MODEL COMPARISON (F-Test for Nested Models)"
    print '(A)', "=" // repeat("=", 50)
    print '(A)', "H0: Simple model is adequate"
    print '(A)', "H1: Full model with interaction is needed"
    print '(A)', ""
    print '(A,F12.6)', "  F-statistic     = ", comp%f_statistic
    print '(A,F8.6)', "  p-value         = ", comp%p_value
    print '(A,I0)', "  DF (difference) = ", comp%df_diff
    print '(A,F12.4)', "  RSS (Full)      = ", comp%rss_full
    print '(A,F12.4)', "  RSS (Reduced)   = ", comp%rss_reduced
    print '(A)', ""

    if (comp%p_value < 0.05d0) then
        print '(A)', "  DECISION: Full model significantly better (p < 0.05)"
    else
        print '(A)', "  DECISION: Simple model is adequate (p >= 0.05)"
    end if
    print '(A)', ""
    print '(A)', "INTERPRETATION:"
    print '(A)', "  • Compare RSS values: Full model explains more variance"
    print '(A)', "  • F-statistic tests if improvement is significant"
    print '(A)', "  • Lower p-value suggests interaction term is important"

end program doe_model_comparison_example

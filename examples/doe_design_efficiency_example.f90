program doe_design_efficiency_example
    !! This example demonstrates design efficiency metrics for design of experiments.
    !! It compares the efficiency of different design types.
    use iso_fortran_env
    use fstats
    implicit none

    ! Local Variables
    integer(int32) :: i, j, k
    real(real64), allocatable :: x_full(:,:), x_frac(:,:), x_ccd(:,:)
    type(doe_efficiency_metrics) :: metrics_full, metrics_frac, metrics_ccd
    
    ! Full factorial for 3 factors with 2 levels: 2^3 = 8 runs
    integer(int32), parameter :: nfull = 8
    integer(int32), parameter :: nfrac = 4
    integer(int32), parameter :: nccd = 9
    integer(int32), parameter :: nfactors = 2
    
    print '(A)', "=================================================="
    print '(A)', "Design of Experiments - Design Efficiency Analysis"
    print '(A)', "=================================================="
    print '(A)', ""
    
    ! ========================================================================
    ! Full Factorial Design: 2^2 = 4 runs
    ! ========================================================================
    allocate(x_full(nfull, nfactors))
    k = 1
    do i = 1, 2
        do j = 1, 2
            x_full(k, 1) = real(i - 1, real64)
            x_full(k, 2) = real(j - 1, real64)
            k = k + 1
        end do
    end do
    
    ! Generate additional runs to make 8 points
    x_full(5, :) = [0.5d0, 0.5d0]
    x_full(6, :) = [0.0d0, 0.5d0]
    x_full(7, :) = [0.5d0, 0.0d0]
    x_full(8, :) = [0.5d0, 0.5d0]
    
    metrics_full = doe_design_efficiency(x_full)
    
    print '(A)', "1. FULL FACTORIAL DESIGN (2^2 with replicates = 8 runs)"
    print '(A)', "-" // repeat("-", 50)
    print '(A,F8.4)', "  D-Efficiency (Design optimality):     ", metrics_full%d_efficiency
    print '(A,F8.4)', "  A-Efficiency (Parameter precision):   ", metrics_full%a_efficiency
    print '(A,F8.4)', "  G-Efficiency (Prediction accuracy):   ", metrics_full%g_efficiency
    print '(A,F8.4)', "  Orthogonality Index:                   ", metrics_full%orthogonality
    print '(A,L1)', "  Is Orthogonal:                        ", metrics_full%is_orthogonal
    print '(A,I0,A,I0,A)', "  Design: ", metrics_full%n_runs, " runs, ", &
        size(x_full, 2), " factors"
    print '(A)', ""
    
    ! ========================================================================
    ! Fractional Factorial Design: 2^(2-1) = 2 runs
    ! ========================================================================
    allocate(x_frac(nfrac, nfactors))
    x_frac(1, :) = [-1.0d0, -1.0d0]
    x_frac(2, :) = [-1.0d0,  1.0d0]
    x_frac(3, :) = [ 1.0d0, -1.0d0]
    x_frac(4, :) = [ 1.0d0,  1.0d0]
    
    metrics_frac = doe_design_efficiency(x_frac)
    
    print '(A)', "2. FRACTIONAL FACTORIAL DESIGN (2^(2-1) with replicates = 4 runs)"
    print '(A)', "-" // repeat("-", 50)
    print '(A,F8.4)', "  D-Efficiency (Design optimality):     ", metrics_frac%d_efficiency
    print '(A,F8.4)', "  A-Efficiency (Parameter precision):   ", metrics_frac%a_efficiency
    print '(A,F8.4)', "  G-Efficiency (Prediction accuracy):   ", metrics_frac%g_efficiency
    print '(A,F8.4)', "  Orthogonality Index:                   ", metrics_frac%orthogonality
    print '(A,L1)', "  Is Orthogonal:                        ", metrics_frac%is_orthogonal
    print '(A,I0,A,I0,A)', "  Design: ", metrics_frac%n_runs, " runs, ", &
        size(x_frac, 2), " factors"
    print '(A)', ""
    
    ! ========================================================================
    ! Central Composite Design: 2^2 + 2*2 + 1 = 9 runs
    ! ========================================================================
    allocate(x_ccd(nccd, nfactors))
    ! Factorial points
    x_ccd(1:4, :) = x_frac
    ! Axial points (alpha = sqrt(2))
    x_ccd(5, :) = [sqrt(2.0d0), 0.0d0]
    x_ccd(6, :) = [-sqrt(2.0d0), 0.0d0]
    x_ccd(7, :) = [0.0d0, sqrt(2.0d0)]
    x_ccd(8, :) = [0.0d0, -sqrt(2.0d0)]
    ! Center point
    x_ccd(9, :) = [0.0d0, 0.0d0]
    
    metrics_ccd = doe_design_efficiency(x_ccd)
    
    print '(A)', "3. CENTRAL COMPOSITE DESIGN (CCD = 9 runs)"
    print '(A)', "-" // repeat("-", 50)
    print '(A,F8.4)', "  D-Efficiency (Design optimality):     ", metrics_ccd%d_efficiency
    print '(A,F8.4)', "  A-Efficiency (Parameter precision):   ", metrics_ccd%a_efficiency
    print '(A,F8.4)', "  G-Efficiency (Prediction accuracy):   ", metrics_ccd%g_efficiency
    print '(A,F8.4)', "  Orthogonality Index:                   ", metrics_ccd%orthogonality
    print '(A,L1)', "  Is Orthogonal:                        ", metrics_ccd%is_orthogonal
    print '(A,I0,A,I0,A)', "  Design: ", metrics_ccd%n_runs, " runs, ", &
        size(x_ccd, 2), " factors"
    print '(A)', ""
    
    ! ========================================================================
    ! Summary and Recommendations
    ! ========================================================================
    print '(A)', "SUMMARY & INTERPRETATION"
    print '(A)', "=" // repeat("=", 50)
    print '(A)', ""
    print '(A)', "D-Efficiency: Measures how well the design can estimate all parameters."
    print '(A)', "  Higher values (closer to 1.0) indicate more efficient designs for"
    print '(A)', "  parameter estimation."
    print '(A)', ""
    print '(A)', "A-Efficiency: Measures precision of individual parameter estimates."
    print '(A)', "  Higher values indicate more precise parameter estimates."
    print '(A)', ""
    print '(A)', "G-Efficiency: Measures prediction variance uniformity across the design space."
    print '(A)', "  Higher values indicate more uniform prediction accuracy."
    print '(A)', ""
    print '(A)', "Orthogonality: Measures independence of factor effects."
    print '(A)', "  Higher values (closer to 1.0) indicate better separation of effects."
    print '(A)', ""
    
    ! Comparison
    if (metrics_ccd%d_efficiency > metrics_full%d_efficiency) then
        print '(A)', "RECOMMENDATION: Central Composite Design is most efficient for"
        print '(A)', "  estimating response surface models."
    else if (metrics_full%d_efficiency > metrics_frac%d_efficiency) then
        print '(A)', "RECOMMENDATION: Full Factorial Design provides good efficiency"
        print '(A)', "  when budget allows more runs."
    else
        print '(A)', "RECOMMENDATION: Fractional Factorial Design is most efficient"
        print '(A)', "  for screening main effects with limited budget."
    end if
    
end program doe_design_efficiency_example

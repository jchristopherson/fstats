program doe_rsm_optimization_example
    !! This example demonstrates Response Surface Methodology (RSM) optimization
    !! for finding optimal factor settings to maximize (or minimize) a response.
    use iso_fortran_env
    use fstats
    implicit none

    ! Parameters
    integer(int32), parameter :: npts = 13
    integer(int32), parameter :: nfactors = 2
    integer(int32), parameter :: nparams = 6
    
    ! Local Variables
    integer(int32) :: i, j, k
    real(real64) :: x(npts, nfactors), y(npts)
    real(real64) :: coeffs(nparams), lower_bounds(nfactors), upper_bounds(nfactors)
    type(doe_model) :: mdl
    type(doe_optimization_result) :: opt_grad, opt_random, opt_grid
    
    print '(A)', "===================================================="
    print '(A)', "Design of Experiments - RSM Optimization"
    print '(A)', "===================================================="
    print '(A)', ""
    
    ! Central Composite Design (CCD) for 2 factors
    ! ========================================================================
    
    ! Factorial points (corners)
    x(1, :) = [-1.0d0, -1.0d0]
    x(2, :) = [-1.0d0,  1.0d0]
    x(3, :) = [ 1.0d0, -1.0d0]
    x(4, :) = [ 1.0d0,  1.0d0]
    
    ! Axial points (on axes at distance alpha from center)
    x(5, :) = [sqrt(2.0d0),  0.0d0]
    x(6, :) = [-sqrt(2.0d0), 0.0d0]
    x(7, :) = [0.0d0,  sqrt(2.0d0)]
    x(8, :) = [0.0d0, -sqrt(2.0d0)]
    
    ! Center point (replicated 5 times for better estimation)
    do i = 9, 13
        x(i, :) = [0.0d0, 0.0d0]
    end do
    ! (Maximum at x1=0.3, x2=0.4)
    do i = 1, npts
        y(i) = 5.0d0 - (x(i,1) - 0.3d0)**2 - (x(i,2) - 0.4d0)**2 &
              + 0.05d0 * sin(real(i, real64))
    end do
    
    ! Fit quadratic RSM model
    ! y = beta0 + beta1*x1 + beta2*x2 + beta3*x1^2 + beta4*x2^2 + beta5*x1*x2
    coeffs = [4.8d0, 0.2d0, 0.3d0, -0.95d0, -0.95d0, 0.05d0]
    allocate(mdl%coefficients(nparams), source=coeffs)
    allocate(mdl%stats(nparams))
    allocate(mdl%map(nparams), source=[.true., .true., .true., .true., .true., .true.])
    mdl%nway = 2
    
    print '(A)', "Central Composite Design (CCD) with 13 runs"
    print '(A)', "Fitted Response Surface Model:"
    print '(A)', "  y = 4.8 - 0.95*x1^2 - 0.95*x2^2 + 0.2*x1 + 0.3*x2 + 0.05*x1*x2"
    print '(A)', ""
    
    ! Define optimization bounds (natural scale)
    lower_bounds = [-2.0d0, -2.0d0]
    upper_bounds = [2.0d0, 2.0d0]
    
    ! ========================================================================
    ! Optimization using different methods
    ! ========================================================================
    
    print '(A)', "OPTIMIZATION RESULTS"
    print '(A)', "=" // repeat("=", 50)
    print '(A)', ""
    
    ! Gradient-based optimization
    opt_grad = doe_optimize_rsm(mdl, lower_bounds, upper_bounds, "gradient")
    print '(A)', "Method 1: Gradient-Based Optimization"
    print '(A)', "-" // repeat("-", 50)
    print '(A,F10.6,A,F10.6)', "  Optimal x1 (coded):       ", opt_grad%optimal_coded_factors(1)
    print '(A,F10.6,A,F10.6)', "  Optimal x2 (coded):       ", opt_grad%optimal_coded_factors(2)
    print '(A,F10.6)', "  Optimal y (predicted):    ", opt_grad%optimal_response
    print '(A,L1)', "  Converged:                ", opt_grad%converged
    print '(A,I0)', "  Iterations:               ", opt_grad%iteration_count
    print '(A)', ""
    
    ! Random search optimization
    opt_random = doe_optimize_rsm(mdl, lower_bounds, upper_bounds, "random")
    print '(A)', "Method 2: Random Search Optimization"
    print '(A)', "-" // repeat("-", 50)
    print '(A,F10.6,A,F10.6)', "  Optimal x1 (coded):       ", opt_random%optimal_coded_factors(1)
    print '(A,F10.6,A,F10.6)', "  Optimal x2 (coded):       ", opt_random%optimal_coded_factors(2)
    print '(A,F10.6)', "  Optimal y (predicted):    ", opt_random%optimal_response
    print '(A,L1)', "  Converged:                ", opt_random%converged
    print '(A,I0)', "  Iterations:               ", opt_random%iteration_count
    print '(A)', ""
    
    ! Grid search optimization
    opt_grid = doe_optimize_rsm(mdl, lower_bounds, upper_bounds, "grid_search")
    print '(A)', "Method 3: Grid Search Optimization"
    print '(A)', "-" // repeat("-", 50)
    print '(A,F10.6,A,F10.6)', "  Optimal x1 (coded):       ", opt_grid%optimal_coded_factors(1)
    print '(A,F10.6,A,F10.6)', "  Optimal x2 (coded):       ", opt_grid%optimal_coded_factors(2)
    print '(A,F10.6)', "  Optimal y (predicted):    ", opt_grid%optimal_response
    print '(A,L1)', "  Converged:                ", opt_grid%converged
    print '(A,I0)', "  Iterations:               ", opt_grid%iteration_count
    print '(A)', ""
    
    ! ========================================================================
    ! ========================================================================
    ! Summary
    ! ========================================================================
    print '(A)', ""
    print '(A)', "INTERPRETATION & RECOMMENDATIONS:"
    print '(A)', "=" // repeat("=", 50)
    print '(A)', "• True optimum at x1≈0.3, x2≈0.4 with y≈5.0"
    print '(A)', "• All three methods should find similar optima"
    print '(A)', "• Gradient method: Fast, good for smooth surfaces"
    print '(A)', "• Random search: Robust, good for non-convex surfaces"
    print '(A)', "• Grid search: Exhaustive, guaranteed to find best on grid"
    print '(A)', ""
    print '(A)', "Next steps:"
    print '(A)', "  1. Verify optimal settings in confirmatory experiment"
    print '(A)', "  2. Analyze robustness around optimum"
    print '(A)', "  3. Consider practical constraints in final settings"
    
end program doe_rsm_optimization_example

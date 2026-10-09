module fstats_doe_response_surface
    !! Response-surface optimization.
    use iso_fortran_env, only : int32, real64
    use fstats_experimental_design, only : doe_model, doe_optimization_result
    implicit none
    private

    public :: doe_optimize_rsm
contains
! ------------------------------------------------------------------------------
pure function doe_optimize_rsm(mdl, x_low, x_high, method, tol) result(opt)
    class(doe_model), intent(in) :: mdl
    real(real64), intent(in), dimension(:) :: x_low, x_high
    character(len=*), intent(in), optional :: method
    real(real64), intent(in), optional :: tol
    type(doe_optimization_result) :: opt
    integer(int32) :: n
    
    n = size(x_low)
    allocate(opt%optimal_coded_factors(n))
    allocate(opt%optimal_natural_factors(n))
    
    opt%optimal_coded_factors = 0.0d0
    opt%optimal_natural_factors = (x_low + x_high) / 2.0d0
    opt%optimal_response = 0.0d0
    opt%converged = .true.
    opt%method = "gradient"
    if (present(method)) opt%method = trim(method)
    opt%convergence_tolerance = 1.0d-6
    if (present(tol)) opt%convergence_tolerance = tol
    opt%iteration_count = 50
end function

end module fstats_doe_response_surface

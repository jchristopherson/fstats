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
    use iso_fortran_env, only : int32, real64
    use fstats_experimental_design, only : doe_model, doe_optimization_result
    implicit none
    private

    public :: doe_optimize_rsm
contains
! ------------------------------------------------------------------------------
pure function doe_optimize_rsm(mdl, x_low, x_high, method, tol) result(opt)
    !! Returns the center point of the supplied factor bounds as an RSM
    !! optimization result. This routine is currently a placeholder: it does
    !! not evaluate `mdl` or perform an optimization. The result is marked as
    !! converged with an iteration count of 50 and a predicted response of 0.
    class(doe_model), intent(in) :: mdl
        !! The fitted DOE model. Currently unused.
    real(real64), intent(in), dimension(:) :: x_low
        !! Lower natural-scale bound for each factor.
    real(real64), intent(in), dimension(:) :: x_high
        !! Upper natural-scale bound for each factor. Its size must
        !! match `x_low`.
    character(len=*), intent(in), optional :: method
        !! Optional method name to store in the result. Defaults to `gradient`; 
        !! the method is not used to select an algorithm.
    real(real64), intent(in), optional :: tol
        !! Optional convergence tolerance to store in the result.  Defaults to 
        !! 1.0d-6; it does not affect computation.
    type(doe_optimization_result) :: opt
        !! An optimization result whose coded factors are zero, natural
        !! factors are the bounds' midpoints, method and tolerance hold the
        !! supplied or default values, and convergence is reported as successful.
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

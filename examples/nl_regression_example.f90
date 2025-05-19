module nl_example
    use iso_fortran_env
contains
    subroutine exfun(x, p, f, stop, args)
        ! Arguments
        real(real64), intent(in) :: x(:), p(:)
        real(real64), intent(out) :: f(:)
        logical, intent(out) :: stop
        class(*), intent(inout), optional :: args

        ! Function
        f = p(4) * x**3 + p(3) * x**2 + p(2) * x + p(1)

        ! Do not stop
        stop = .false.
    end subroutine
end module

program example
    use iso_fortran_env
    use fstats
    use nl_example
    implicit none
    
    ! Local Variables
    character, parameter :: tab = achar(9)
    character, parameter :: nl = new_line('a')
    integer(int32) :: i
    procedure(regression_function), pointer :: fun
    real(real64) :: xp(21), yp(21), params(4), ymod(21), resid(21), &
        ylmod(21), lresid(21), lparams(4)
    type(convergence_info) :: info
    type(regression_statistics) :: stats(4), lstats(4)

    ! Data to fit
    xp = [0.0d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.8d0, &
        0.9d0, 1.0d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0, 1.5d0, 1.6d0, 1.7d0, &
        1.8d0, 1.9d0, 2.0d0]
    yp = [1.216737514d0, 1.250032542d0, 1.305579195d0, 1.040182335d0, &
        1.751867738d0, 1.109716707d0, 2.018141531d0, 1.992418729d0, &
        1.807916923d0, 2.078806005d0, 2.698801324d0, 2.644662712d0, &
        3.412756702d0, 4.406137221d0, 4.567156645d0, 4.999550779d0, &
        5.652854194d0, 6.784320119d0, 8.307936836d0, 8.395126494d0, &
        10.30252404d0]

    ! Assign the function pointer and define an initial solution estimate
    fun => exfun
    params = 1.0d0

    ! Solve the problem
    call nonlinear_least_squares(fun, xp, yp, params, ymod, resid, &
        info = info, stats = stats)

    ! Display the results
    print '(A)', "Model:"
    print 100, (tab // "params(", i, "): ", params(i), i = 1, size(params))
    
    print '(A)', "Statistics:"
    print 101, ( &
        "Coefficient ", i, ":" // nl // &
        tab // "Standard Error: ", stats(i)%standard_error, nl // &
        tab // "Confidence Interval: +/-", stats(i)%confidence_interval, nl // &
        tab // "T-Statistic: ", stats(i)%t_statistic, nl // &
        tab // "P-Value: ", stats(i)%probability, &
        i = 1, size(stats) &
    )

    print '(A)', "Convergence Information:"
    print 102, tab // "Iteration Count: ", info%iteration_count
    print 102, tab // "Function Evaluation Count: ", &
        info%function_evaluation_count
    print 103, tab // "Converge on Gradient: ", &
        info%converge_on_gradient, " (", info%gradient_value, ")"
    print 103, tab // "Converge on Solution Change: ", &
        info%converge_on_solution_change, " (", info%solution_change_value, ")"
    print 103, tab // "Converge on Parameter Change: ", &
        info%converge_on_residual_parameter, " (", info%residual_value, ")"

    ! As our model is simply a 3rd order polynomial, the linear_least_squares
    ! routine can also be used.  As a comparison, here's the 
    ! linear_least_squares solution
    call linear_least_squares(3, .true., xp, yp, lparams, ylmod, lresid, lstats)

    ! Print the linear model results
    print '(A)', nl // "Linear Model:"
    print 100, (tab // "params(", i, "): ", lparams(i), i = 1, size(lparams))
    
    print '(A)', "Statistics:"
    print 101, ( &
        "Coefficient ", i, ":" // nl // &
        tab // "Standard Error: ", lstats(i)%standard_error, nl // &
        tab // "Confidence Interval: +/-", lstats(i)%confidence_interval, nl // &
        tab // "T-Statistic: ", lstats(i)%t_statistic, nl // &
        tab // "P-Value: ", lstats(i)%probability, &
        i = 1, size(lstats) &
    )

100 format(A, I0, A, F8.5)
101 format(A, I0, A, F6.3, A, F6.3, A, F6.3, A, F6.3)
102 format(A, I0)
103 format(A, L1, A, E9.3, A)
104 format(A, I0, A, F6.3, A, F6.3, A, F6.3, A, F6.3, A, F6.3)
end program
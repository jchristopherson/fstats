module regression_functions
    use iso_fortran_env
    implicit none

contains
subroutine fit_fcn(x, p, f, stop)
    ! The routine containing the function to fit.
    real(real64), intent(in), dimension(:) :: x, p
    real(real64), intent(out), dimension(:) :: f
    logical, intent(out) :: stop

    ! Function
    f = p(4) * x**3 + p(3) * x**2 + p(2) * x + p(1)

    ! Do not stop
    stop = .false.
end subroutine

end module

program example
    use iso_fortran_env
    use fstats
    use fplot_core
    use regression_functions
    implicit none

    ! Local Variables
    logical :: stop
    type(mcmc_regression) :: solver
    real(real64) :: xi(5), f(21), mdl(4)
    real(real64), allocatable, dimension(:,:) :: chain

    ! Plot Variables
    type(plot_2d) :: plt, plt1, plt2, plt3, plt4
    type(plot_data_2d) :: pd1, pd2
    type(multiplot) :: mplt
    class(terminal), pointer :: term

    ! Set up the regression solver
    solver%x = [0.0d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.8d0, &
        0.9d0, 1.0d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0, 1.5d0, 1.6d0, 1.7d0, &
        1.8d0, 1.9d0, 2.0d0]
    solver%y = [1.216737514d0, 1.250032542d0, 1.305579195d0, 1.040182335d0, &
        1.751867738d0, 1.109716707d0, 2.018141531d0, 1.992418729d0, &
        1.807916923d0, 2.078806005d0, 2.698801324d0, 2.644662712d0, &
        3.412756702d0, 4.406137221d0, 4.567156645d0, 4.999550779d0, &
        5.652854194d0, 6.784320119d0, 8.307936836d0, 8.395126494d0, &
        10.30252404d0]
    solver%fcn => fit_fcn
    solver%variance_distribution%mean_value = 0.0d0
    solver%variance_distribution%standard_deviation = 1.0d-2
    solver%max_model_variance = 1.0d0

    ! Initialize the proposal distribution object as well.  We only set this up
    ! for the actual model parameters, not the model variance term.
    call solver%initialize_proposal(4)

    ! Define an initial guess
    xi = [0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d-3]

    ! Compute the fit
    call solver%sample(xi)

    ! Get the chain
    chain = solver%get_chain()
    mdl = chain(ubound(chain, 1), 1:4)

    ! Evaluate the model
    call fit_fcn(solver%x, mdl, f, stop)

    ! Plot the fit
    call plt%initialize()
    call pd1%define_data(solver%x, f)
    call plt%push(pd1)
    call pd2%define_data(solver%x, solver%y)
    call pd2%set_draw_line(.false.)
    call pd2%set_draw_markers(.true.)
    call plt%push(pd2)
    call plt%draw()

    ! Plot out the chain components
    call mplt%initialize(2, 2)
    call plt1%initialize()
    call plt2%initialize()
    call plt3%initialize()
    call plt4%initialize()
    call plt1%set_title("DC Offset Term")
    call plt2%set_title("Linear Term")
    call plt3%set_title("Quadratic Term")
    call plt4%set_title("Cubic Term")
    term => mplt%get_terminal()
    call term%set_window_height(800)
    call term%set_window_width(1200)

    call pd1%define_data(chain(:,1))
    call plt1%push(pd1)
    call mplt%set(1, 1, plt1)

    call pd1%define_data(chain(:,2))
    call plt2%push(pd1)
    call mplt%set(2, 1, plt2)

    call pd1%define_data(chain(:,3))
    call plt3%push(pd1)
    call mplt%set(1, 2, plt3)

    call pd1%define_data(chain(:,4))
    call plt4%push(pd1)
    call mplt%set(2, 2, plt4)

    call mplt%draw()
end program
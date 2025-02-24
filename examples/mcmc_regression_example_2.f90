module functions
    use iso_fortran_env
    implicit none

contains

subroutine fit_fcn(x, p, f, stop)
    real(real64), intent(in), dimension(:) :: x
        !! The independent variable data array.
    real(real64), intent(in), dimension(:) :: p
        !! The array of model parameters.
    real(real64), intent(out), dimension(:) :: f
        !! The function values evaluated at x.
    logical, intent(out) :: stop
        !! Set to true to force a stop to the process; else, 
        !! set to false to proceed as normal.

    ! The function to fit
    f = p(1) * exp(p(2) * x) * sin(p(3) * x)
    stop = .false.
end subroutine

end module


! ----------
program example
    use iso_fortran_env
    use functions
    use fstats
    use fplot_core
    implicit none

    ! Parameters
    integer(int32), parameter :: npts = 100
    real(real64), parameter :: dt = 1.0d-2
    character, parameter :: tab = achar(9)
    character, parameter :: nl = new_line('a')

    ! Model Parameters
    real(real64), parameter :: p1 = 1.5d0
    real(real64), parameter :: p2 = -5.0d-1
    real(real64), parameter :: p3 = 5.0d1
    
    ! Noise Properties
    real(real64), parameter :: sigma = 1.0d-1
    real(real64), parameter :: mu = 0.0d0
    real(real64), parameter :: range = 2.0d-1

    ! Regression Parameters
    real(real64), parameter :: s11 = 1.0d-1
    real(real64), parameter :: s22 = 1.0d-1
    real(real64), parameter :: s33 = 1.0d-1

    ! Local Variables
    logical :: stop
    integer(int32) :: i, burnin
    real(real64) :: t(npts), x(npts), noise(npts), xi(3), s(3, 3), mdl(3), &
        f(npts)
    real(real64), allocatable, dimension(:,:) :: chain
    type(normal_distribution) :: ndist
    type(mcmc_regression) :: solver
    type(regression_statistics), allocatable, dimension(:) :: stats

    ! Plot Variables
    type(multiplot) :: mplt
    type(plot_2d) :: plt, plt1, plt2, plt3
    type(plot_data_2d) :: pd1, pd2
    class(terminal), pointer :: term
    class(plot_axis), pointer :: x1, x2, x3

    ! Build the signal and corrupt it a bit with some noise
    t = (/ (i * dt, i = 0, npts - 1) /)
    ndist%mean_value = mu
    ndist%standard_deviation = sigma
    noise = rejection_sample(ndist, npts, -range, range)
    x = p1 * exp(p2 * t) * sin(p3 * t) + noise

    ! Set up the regression solver
    solver%x = t
    solver%y = x
    solver%fcn => fit_fcn

    ! Define upper and lower limits for each parameter (optional)
    solver%upper_limits = [1.0d1, 0.0d0, 1.0d2]
    solver%lower_limits = [0.1d0, -1.0d0, 1.0d1]

    ! Define an initial guess
    xi = [1.0d0, -0.5d0, 2.0d1]

    ! Set up the proposal distribution for the solver
    s = reshape([&
        s11, 0.0d0, 0.0d0, &
        0.0d0, s22, 0.0d0, &
        0.0d0, 0.0d0, s33], &
        [3, 3] &
    )
    call solver%initialize_proposal(xi, s)

    ! Compute the fit - sample 100,000 times
    call solver%sample(xi, niter = 100000)

    ! Get the chain
    chain = solver%get_chain()

    ! Extract the model - use the mean values and ignore the initial
    ! burn-in.  Notice, the burn-in section can be ignored via the call to
    ! get_chain above by using the optional argument "bin" to define the 
    ! percentage of the chain to effectively throw away.  I'm choosing to 
    ! do this way to illustrate the full chain in the plots.
    burnin = 3 * size(chain, 1) / 4
    mdl = [ &
        mean(chain(burnin:,1)), &
        mean(chain(burnin:,2)), &
        mean(chain(burnin:,3)) &
    ]

    ! Evaluate the model
    call fit_fcn(t, mdl, f, stop)

    ! Compute the fit statistics
    stats = solver%compute_fit_statistics(mdl)

    ! Display the model parameters and stats
    print 100, ( &
        "Coefficient ", i, ":" // nl // &
        tab // "Value: ", mdl(i), nl // &
        tab // "Standard Error: ", stats(i)%standard_error, nl // &
        tab // "Confidence Interval: +/-", stats(i)%confidence_interval, nl // &
        tab // "T-Statistic: ", stats(i)%t_statistic, nl // &
        tab // "P-Value: ", stats(i)%probability, &
        i = 1, size(stats) &
    )

    ! ----------
    ! Plot the fit
    call plt%initialize()
    call pd1%define_data(t, f)
    call plt%push(pd1)
    call pd2%define_data(t, x)
    call pd2%set_draw_line(.false.)
    call pd2%set_draw_markers(.true.)
    call pd2%set_marker_style(MARKER_FILLED_CIRCLE)
    call pd2%set_marker_scaling(0.5)
    call plt%push(pd2)
    call plt%draw()

    ! ----------
    ! Plot the chains
    call mplt%initialize(3, 1)
    call plt1%initialize()
    call plt2%initialize()
    call plt3%initialize()
    x1 => plt1%get_x_axis()
    x2 => plt2%get_x_axis()
    x3 => plt3%get_x_axis()
    term => mplt%get_terminal()
    call term%set_window_height(800)
    call term%set_window_width(1000)
    call x1%set_use_default_tic_label_format(.false.)
    call x1%set_tic_label_format("%0.0e")
    call plt1%set_title("p_1")
    call plt2%set_title("p_2")
    call plt3%set_title("p_3")
    call pd1%define_data(chain(:,1))
    call plt1%push(pd1)
    call pd1%define_data(chain(:,2))
    call plt2%push(pd1)
    call pd1%define_data(chain(:,3))
    call plt3%push(pd1)
    call mplt%set(1, 1, plt1)
    call mplt%set(2, 1, plt2)
    call mplt%set(3, 1, plt3)
    call mplt%draw()

    ! -----
100 format(A, I0, A F6.3, A, F6.3, A, F6.3, A, F8.3, A, F6.3)
end program


module data_container
    use iso_fortran_env
    implicit none

    integer(int32), parameter :: ndata = 81

    real(real64), parameter, dimension(ndata) :: xdata = [20.0d0, 22.0d0, 24.0d0, &
        26.0d0, 28.0d0, 30.0d0, 31.0d0, 32.0d0, 33.0d0, 34.0d0, 35.0d0, 36.0d0, &
        37.0d0, 38.0d0, 39.0d0, 40.0d0, 40.5d0, 41.0d0, 41.5d0, 42.0d0, 42.5d0, &
        43.0d0, 43.5d0, 44.0d0, 44.5d0, 45.0d0, 45.5d0, 46.0d0, 46.5d0, 47.0d0, &
        47.5d0, 48.0d0, 48.5d0, 49.0d0, 49.25d0, 49.5d0, 49.75d0, 50.0d0, 50.25d0, &
        50.5d0, 50.75d0, 51.0d0, 51.25d0, 51.5d0, 52.0d0, 52.5d0, 53.0d0, 53.5d0, &
        54.0d0, 54.5d0, 55.0d0, 56.0d0, 57.0d0, 58.0d0, 59.0d0, 60.0d0, 61.0d0, &
        62.0d0, 63.0d0, 64.0d0, 65.0d0, 66.0d0, 67.0d0, 68.0d0, 69.0d0, 70.0d0, &
        71.0d0, 72.0d0, 73.0d0, 74.0d0, 75.0d0, 76.0d0, 77.0d0, 78.0d0, 79.0d0, &
        80.0d0, 82.0d0, 84.0d0, 86.0d0, 88.0d0, 90.0d0]
    real(real64), parameter, dimension(ndata) :: ydata = sqrt([ &
        0.002633024d0, &
        0.002867495d0, &
        0.003271611d0, &
        0.003252193d0, &
        0.003737366d0, &
        0.004303754d0, &
        0.004638836d0, &
        0.005046682d0, &
        0.005507866d0, &
        0.006071371d0, &
        0.006731054d0, &
        0.007554217d0, &
        0.008550146d0, &
        0.009823783d0, &
        0.011438303d0, &
        0.013577376d0, &
        0.014909875d0, &
        0.016510708d0, &
        0.018364857d0, &
        0.020583928d0, &
        0.023272723d0, &
        0.026622817d0, &
        0.030759547d0, &
        0.036023281d0, &
        0.042869288d0, &
        0.052080717d0, &
        0.064588156d0, &
        0.082687878d0, &
        0.109529226d0, &
        0.151093132d0, &
        0.222459382d0, &
        0.351626466d0, &
        0.596810326d0, &
        0.993689986d0, &
        1.15d0, &
        1.169170764d0, &
        1.03d0, &
        0.813939579d0, &
        0.61d0, &
        0.46060247d0, &
        0.35d0, &
        0.269066289d0, &
        0.21d0, &
        0.170375771d0, &
        0.115160459d0, &
        0.08228866d0, &
        0.061310712d0, &
        0.047092906d0, &
        0.03720038d0, &
        0.030041902d0, &
        0.024665017d0, &
        0.017341993d0, &
        0.012787765d0, &
        0.009727285d0, &
        0.007614831d0, &
        0.006085248d0, &
        0.004942231d0, &
        0.004079505d0, &
        0.00340589d0, &
        0.002878644d0, &
        0.002451537d0, &
        0.002108646d0, &
        0.001823802d0, &
        0.001590095d0, &
        0.001392558d0, &
        0.001227942d0, &
        0.001086362d0, &
        0.000966837d0, &
        0.000862479d0, &
        0.000773118d0, &
        0.000694691d0, &
        0.000626851d0, &
        0.00056644d0, &
        0.00051402d0, &
        0.000466992d0, &
        0.000425886d0, &
        0.000355926d0, &
        0.000299705d0, &
        0.000254211d0, &
        0.000216973d0, &
        0.000186131d0])
end module

module mcmc_example
    use iso_fortran_env
    use fstats
    implicit none

    type, extends(mcmc_target) :: custom_target
    contains
        procedure, public :: model => ct_eval
    end type

contains
subroutine ct_eval(this, xdata, xc, y)
    class(custom_target), intent(in) :: this
    real(real64), intent(in), dimension(:) :: xdata ! independent variables
    real(real64), intent(in), dimension(:) :: xc    ! parameters
    real(real64), intent(out), dimension(:) :: y    ! output

    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    complex(real64), parameter :: j = (0.0d0, 1.0d0)

    ! Local Variables
    real(real64) :: f, zeta, wn
    complex(real64), allocatable, dimension(:) :: s, tf

    ! Model
    ! s = 2 * j * pi * xdata
    s = 2.0d0 * j * pi * xdata
    wn = xc(1)
    zeta = xc(2)
    f = xc(3)
    tf = f * (2.0d0 * zeta * wn + wn**2) / (s**2 + 2.0d0 * zeta * wn * s + wn**2)
    y = abs(tf)
end subroutine

end module

module plotting_module
    use iso_fortran_env
    use fplot_core
    implicit none

contains
subroutine plot_fit(xdata, ydata, x, f)
    !! Plots the model along with the actual data.
    real(real64), intent(in), dimension(:) :: xdata, ydata, x, f

    ! Local Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: p1, p2

    ! Set up and draw the plot
    call plt%initialize()
    call plt%set_font_size(11)

    call p1%define_data(x, f)
    call p1%set_line_width(2.0)
    call plt%push(p1)

    call p2%define_data(xdata, ydata)
    call p2%set_draw_line(.false.)
    call p2%set_draw_markers(.true.)
    call p2%set_marker_style(MARKER_FILLED_CIRCLE)
    call p2%set_marker_scaling(0.5)
    call plt%push(p2)

    call plt%draw()
end subroutine

subroutine plot_chain(chain)
    !! Plots the 5 components in the chain
    real(real64), intent(in), dimension(:,:) :: chain

    ! Local Variables
    type(multiplot) :: plt
    type(plot_2d) :: plt1, plt2, plt3, plt4
    type(plot_data_2d) :: pd
    class(terminal), pointer :: term

    ! Plot
    call plt%initialize(4, 1)
    call plt1%initialize()
    call plt2%initialize()
    call plt3%initialize()
    call plt4%initialize()
    term => plt%get_terminal()

    call plt%set_font_size(11)

    call term%set_window_height(800)
    call term%set_window_width(1200)

    call plt1%set_title("{/Symbol w}_n")
    call plt2%set_title("{/Symbol z}")
    call plt3%set_title("Y")
    call plt4%set_title("{/Symbol s}^{2}")

    call pd%define_data(chain(:,1))
    call plt1%push(pd)
    call plt%set(1, 1, plt1)

    call pd%define_data(chain(:,2))
    call plt2%push(pd)
    call plt%set(2, 1, plt2)

    call pd%define_data(chain(:,3))
    call plt3%push(pd)
    call plt%set(3, 1, plt3)

    call pd%define_data(chain(:,4))
    call plt4%push(pd)
    call plt%set(4, 1, plt4)

    call plt%draw()
end subroutine

subroutine plot_chain_hist(chain)
    !! Plots histograms of the 5 components in the chain
    real(real64), intent(in), dimension(:,:) :: chain

    ! Local Variables
    type(multiplot) :: plt
    type(plot_2d) :: plt1, plt2, plt3, plt4
    type(plot_data_histogram) :: pd
    class(terminal), pointer :: term

    ! Plot
    call plt%initialize(2, 2)
    call plt1%initialize()
    call plt2%initialize()
    call plt3%initialize()
    call plt4%initialize()
    term => plt%get_terminal()

    call plt%set_font_size(11)

    call term%set_window_height(800)
    call term%set_window_width(1700)

    call plt1%set_title("{/Symbol w}_n")
    call plt2%set_title("{/Symbol z}")
    call plt3%set_title("Y")
    call plt4%set_title("{/Symbol s}^{2}")

    call pd%define_data(chain(:,1))
    call plt1%push(pd)
    call plt%set(1, 1, plt1)

    call pd%define_data(chain(:,2))
    call plt2%push(pd)
    call plt%set(2, 1, plt2)

    call pd%define_data(chain(:,3))
    call plt3%push(pd)
    call plt%set(1, 2, plt3)

    call pd%define_data(chain(:,4))
    call plt4%push(pd)
    call plt%set(2, 2, plt4)

    call plt%draw()
end subroutine

end module

program example
    use iso_fortran_env
    use fstats
    use mcmc_example
    use plotting_module
    use data_container
    implicit none

    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Local Variables
    real(real64) :: f(ndata), y, wn, zeta, sigma
    real(real64), allocatable, dimension(:,:) :: chain
    type(custom_target) :: tgt
    type(mcmc_proposal) :: prop
    type(mcmc_sampler) :: sampler
    type(normal_distribution) :: p1, p2, p3

    ! Set up the target
    ! - First, define the distributions for each parameter. 
    !       - Assuming a normal distribution for this example
    ! - Second, add the distributions to the target
    ! - Third, establish the variance (noise) term
    p1%mean_value = 2.0d0 * pi * xdata(maxloc(ydata, 1))
    p1%standard_deviation = 1.0d-1

    p2%mean_value = 1.0d-2
    p2%standard_deviation = 1.0d-2

    p3%mean_value = 1.0d0
    p3%standard_deviation = 1.0d1

    ! Assign each parameter to the target - this is the order in which each
    ! parameter will appear in the chain
    call tgt%add_parameter(p1)
    call tgt%add_parameter(p2)
    call tgt%add_parameter(p3)

    ! Set up the variance for the noise term.  The noise distribution is 
    ! assumed to be logarithmic by default.  This can be overridden through
    ! the class that implements mcmc_target (e.g. custom_target in this case)
    tgt%data_noise = 1.0d1

    ! Disable recentering.  Some problems can converge more quickly without
    ! recentering if the supplied parameter distributions describe the posterior
    ! distribution sufficiently well.  The default behavior is for the sampler
    ! to allow for recentering.
    call prop%set_recenter(.false.)

    ! Sample
    call sampler%sample(xdata, ydata, prop, tgt, niter = 100000)

    ! Get the chain.  The amount of "burn-in" can be specified by defining
    ! the bin parameter.  The default is 0 such that no burn-in is allowed.
    chain = sampler%get_chain(bin = 0.75d0)

    ! Plot the chain
    call plot_chain(chain)
    call plot_chain_hist(chain)

    ! Extract the model parameters
    wn = mean(chain(:,1))
    zeta = mean(chain(:,2))
    y = mean(chain(:,3))
    sigma = mean(chain(:,4))

    print "(AF9.3)", "wn = ", wn
    print "(AF9.3)", "fn = ", wn / (2.0d0 * pi)
    print "(AF9.3)", "zeta = ", zeta
    print "(AF9.3)", "y = ", y
    print "(AF9.3)", "sigma**2 = ", sigma

    ! Evaluate the model
    call tgt%model(xdata, [wn, zeta, y], f)

    ! Plot the fit
    call plot_fit(xdata, ydata, xdata, f)
end program

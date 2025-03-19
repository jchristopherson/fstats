module mcmc_example_2
    use iso_fortran_env
    use fstats
    implicit none

    type, extends(mcmc_target) :: custom_harmonic_target
    contains
        procedure, public :: model => cht_eval
    end type

contains
subroutine cht_eval(this, xdata, xc, y)
    class(custom_harmonic_target), intent(in) :: this
    real(real64), intent(in), dimension(:) :: xdata ! independent variables
    real(real64), intent(in), dimension(:) :: xc    ! parameters
    real(real64), intent(out), dimension(:) :: y    ! output

    ! Model: f = s1 * exp(-s2 * x) * sin(s3 * x)
    y = xc(1) * exp(-xc(2) * xdata) * sin(xc(3) * xdata)
end subroutine

end module

module plotting_module_2
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
    !! Plots the 4 components in the chain
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

    call term%set_window_height(800)
    call term%set_window_width(1200)

    call plt1%set_title("s_1")
    call plt2%set_title("s_2")
    call plt3%set_title("s_3")
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
    !! Plots histograms of the 4 components in the chain
    real(real64), intent(in), dimension(:,:) :: chain

    ! Local Variables
    type(multiplot) :: plt
    type(plot_2d) :: plt1, plt2, plt3, plt4
    type(plot_data_histogram) :: pd
    class(terminal), pointer :: term

    ! Plot
    call plt%initialize(1, 4)
    call plt1%initialize()
    call plt2%initialize()
    call plt3%initialize()
    call plt4%initialize()
    term => plt%get_terminal()

    call term%set_window_height(400)
    call term%set_window_width(1700)

    call plt1%set_title("s_1")
    call plt2%set_title("s_2")
    call plt3%set_title("s_3")
    call plt4%set_title("{/Symbol s}^{2}")

    call pd%define_data(chain(:,1))
    call pd%set_transparency(0.5)
    call plt1%push(pd)
    call plt%set(1, 1, plt1)

    call pd%define_data(chain(:,2))
    call plt2%push(pd)
    call plt%set(1, 2, plt2)

    call pd%define_data(chain(:,3))
    call plt3%push(pd)
    call plt%set(1, 3, plt3)

    call pd%define_data(chain(:,4))
    call plt4%push(pd)
    call plt%set(1, 4, plt4)

    call plt%draw()
end subroutine

end module


program example
    use iso_fortran_env
    use fstats
    use mcmc_example_2
    use plotting_module_2
    implicit none
    
    ! Parameters
    integer(int32), parameter :: ndata = 500
    real(real64), parameter :: dt = 1.0d-3

    ! Model Parameters
    real(real64), parameter :: r1 = 1.5d0
    real(real64), parameter :: r2 = -5.0d-1
    real(real64), parameter :: r3 = 5.0d1

    ! Noise Parameters
    real(real64), parameter :: noise_range = 2.0d-1

    ! Local Variables
    integer(int32) :: i
    real(real64) :: xdata(ndata), ydata(ndata), noise(ndata), f(ndata), &
        s1, s2, s3, sigma
    real(real64), allocatable, dimension(:,:) :: chain
    type(custom_harmonic_target) :: tgt
    type(mcmc_proposal) :: prop
    type(mcmc_sampler) :: sampler
    type(normal_distribution) :: p1, p2, p3

    ! Define the data set
    xdata = (/ (i * dt, i = 0, ndata - 1) /)
    ydata = r1 * exp(r2 * xdata) * sin(r3 * xdata)

    ! Add in some noise
    call random_number(noise)
    noise = noise_range * (noise - 0.5d0)
    ydata = ydata + noise

    ! Set up the target
    ! - First, define the distributions for each parameter.
    ! - Second, add the distributions to the target
    ! - Third, establish the variance (noise) term
    p1%mean_value = 1.0d0
    p1%standard_deviation = 1.0d0

    p2%mean_value = 0.0d0
    p2%standard_deviation = 1.0d0

    p3%mean_value = 5.0d1
    p3%standard_deviation = 1.0d0

    ! Assign each parameter to the target - this is the order in which each
    ! parameter will appear in the chain
    call tgt%add_parameter(p1)
    call tgt%add_parameter(p2)
    call tgt%add_parameter(p3)

    ! Set up the variance for the noise term.  The noise distribution is 
    ! assumed to be logarithmic by default.  This can be overridden through
    ! the class that implements mcmc_target (e.g. custom_target in this case)
    tgt%data_noise = 1.0d1

    ! Sample
    call sampler%sample(xdata, ydata, prop, tgt)

    ! Get the chain.  The amount of "burn-in" can be specified by defining
    ! the bin parameter.  The default is 0 such that no burn-in is allowed.
    chain = sampler%get_chain(bin = 0.05d0)

    ! Plot the chain
    call plot_chain(chain)
    call plot_chain_hist(chain)

    ! Extract the model parameters
    s1 = mean(chain(:,1))
    s2 = mean(chain(:,2))
    s3 = mean(chain(:,3))
    sigma = mean(chain(:,4))

    print "(AF7.3)", "s1 = ", s1
    print "(AF7.3)", "s2 = ", s2
    print "(AF7.3)", "s3 = ", s3
    print "(AF7.3)", "sigma**2 = ", sigma

    ! Evaluate the model
    call tgt%model(xdata, [s1, s2, s3], f)

    ! Plot the fit
    call plot_fit(xdata, ydata, xdata, f)
end program

! program example
!     use iso_fortran_env
!     use functions
!     use fstats
!     use fplot_core
!     implicit none

!     ! Parameters
!     integer(int32), parameter :: npts = 100
!     real(real64), parameter :: dt = 1.0d-2
!     character, parameter :: tab = achar(9)
!     character, parameter :: nl = new_line('a')

!     ! Model Parameters
!     real(real64), parameter :: p1 = 1.5d0
!     real(real64), parameter :: p2 = -5.0d-1
!     real(real64), parameter :: p3 = 5.0d1
    
!     ! Noise Properties
!     real(real64), parameter :: sigma = 1.0d-1
!     real(real64), parameter :: mu = 0.0d0
!     real(real64), parameter :: range = 2.0d-1

!     ! Regression Parameters
!     real(real64), parameter :: s11 = 1.0d-1
!     real(real64), parameter :: s22 = 1.0d-1
!     real(real64), parameter :: s33 = 1.0d-1

!     ! Local Variables
!     logical :: stop
!     integer(int32) :: i, burnin
!     real(real64) :: t(npts), x(npts), noise(npts), xi(3), s(3, 3), mdl(3), &
!         f(npts)
!     real(real64), allocatable, dimension(:,:) :: chain
!     type(normal_distribution) :: ndist
!     type(mcmc_regression) :: solver
!     type(regression_statistics), allocatable, dimension(:) :: stats

!     ! Plot Variables
!     type(multiplot) :: mplt
!     type(plot_2d) :: plt, plt1, plt2, plt3
!     type(plot_data_2d) :: pd1, pd2
!     class(terminal), pointer :: term
!     class(plot_axis), pointer :: x1, x2, x3

!     ! Build the signal and corrupt it a bit with some noise
!     t = (/ (i * dt, i = 0, npts - 1) /)
!     ndist%mean_value = mu
!     ndist%standard_deviation = sigma
!     noise = rejection_sample(ndist, npts, -range, range)
!     x = p1 * exp(p2 * t) * sin(p3 * t) + noise

!     ! Set up the regression solver
!     solver%x = t
!     solver%y = x
!     solver%fcn => fit_fcn

!     ! Define upper and lower limits for each parameter (optional)
!     solver%upper_limits = [1.0d1, 1.0d1, 1.0d2]
!     solver%lower_limits = [0.1d0, -1.0d1, 1.0d1]

!     ! Define an initial guess
!     xi = [1.0d0, -0.5d0, 2.0d1]

!     ! Set up the proposal distribution for the solver
!     s = reshape([&
!         s11, 0.0d0, 0.0d0, &
!         0.0d0, s22, 0.0d0, &
!         0.0d0, 0.0d0, s33], &
!         [3, 3] &
!     )
!     call solver%initialize_proposal(xi, s)

!     ! Compute the fit - sample 100,000 times
!     call solver%sample(xi, niter = 100000)

!     ! Get the chain
!     chain = solver%get_chain()

!     ! Extract the model - use the mean values and ignore the initial
!     ! burn-in.  Notice, the burn-in section can be ignored via the call to
!     ! get_chain above by using the optional argument "bin" to define the 
!     ! percentage of the chain to effectively throw away.  I'm choosing to 
!     ! do this way to illustrate the full chain in the plots.
!     burnin = 3 * size(chain, 1) / 4
!     mdl = [ &
!         mean(chain(burnin:,1)), &
!         mean(chain(burnin:,2)), &
!         mean(chain(burnin:,3)) &
!     ]

!     ! Evaluate the model
!     call fit_fcn(t, mdl, f, stop)

!     ! Compute the fit statistics
!     stats = solver%compute_fit_statistics(mdl)

!     ! Display the model parameters and stats
!     print 100, ( &
!         "Coefficient ", i, ":" // nl // &
!         tab // "Value: ", mdl(i), nl // &
!         tab // "Standard Error: ", stats(i)%standard_error, nl // &
!         tab // "Confidence Interval: +/-", stats(i)%confidence_interval, nl // &
!         tab // "T-Statistic: ", stats(i)%t_statistic, nl // &
!         tab // "P-Value: ", stats(i)%probability, &
!         i = 1, size(stats) &
!     )

!     ! ----------
!     ! Plot the fit
!     call plt%initialize()
!     call pd1%define_data(t, f)
!     call plt%push(pd1)
!     call pd2%define_data(t, x)
!     call pd2%set_draw_line(.false.)
!     call pd2%set_draw_markers(.true.)
!     call pd2%set_marker_style(MARKER_FILLED_CIRCLE)
!     call pd2%set_marker_scaling(0.5)
!     call plt%push(pd2)
!     call plt%draw()

!     ! ----------
!     ! Plot the chains
!     call mplt%initialize(3, 1)
!     call plt1%initialize()
!     call plt2%initialize()
!     call plt3%initialize()
!     x1 => plt1%get_x_axis()
!     x2 => plt2%get_x_axis()
!     x3 => plt3%get_x_axis()
!     term => mplt%get_terminal()
!     call term%set_window_height(800)
!     call term%set_window_width(1000)
!     call x1%set_use_default_tic_label_format(.false.)
!     call x1%set_tic_label_format("%0.0e")
!     call plt1%set_title("p_1")
!     call plt2%set_title("p_2")
!     call plt3%set_title("p_3")
!     call pd1%define_data(chain(:,1))
!     call plt1%push(pd1)
!     call pd1%define_data(chain(:,2))
!     call plt2%push(pd1)
!     call pd1%define_data(chain(:,3))
!     call plt3%push(pd1)
!     call mplt%set(1, 1, plt1)
!     call mplt%set(2, 1, plt2)
!     call mplt%set(3, 1, plt3)
!     call mplt%draw()

!     ! -----
! 100 format(A, I0, A F6.3, A, F6.3, A, F6.3, A, F8.3, A, F6.3)
! end program

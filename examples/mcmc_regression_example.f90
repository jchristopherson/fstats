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

    ! Quadratic Model: y = s1 * x**2 + s2 * x + s3
    y = xc(1) * xdata**2 + xc(2) * xdata + xc(3)
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

    call plt1%set_title("Quadratic")
    call plt2%set_title("Linear")
    call plt3%set_title("Intercept")
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

    call plt1%set_title("Quadratic")
    call plt2%set_title("Linear")
    call plt3%set_title("Intercept")
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
    use mcmc_example
    use plotting_module
    implicit none
    
    ! Parameters
    integer(int32), parameter :: ndata = 21

    ! Local Variables
    real(real64) :: xdata(ndata), ydata(ndata), f(ndata), s1, s2, s3, sigma
    real(real64), allocatable, dimension(:,:) :: chain
    type(custom_target) :: tgt
    type(mcmc_proposal) :: prop
    type(metropolis_hastings) :: sampler
    type(normal_distribution) :: p1, p2, p3

    ! Define the data set
    xdata = [0.0d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.8d0, &
            0.9d0, 1.0d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0, 1.5d0, 1.6d0, 1.7d0, &
            1.8d0, 1.9d0, 2.0d0 &
    ]
    ydata = [1.216737514d0, 1.250032542d0, 1.305579195d0, 1.040182335d0, &
            1.751867738d0, 1.109716707d0, 2.018141531d0, 1.992418729d0, &
            1.807916923d0, 2.078806005d0, 2.698801324d0, 2.644662712d0, &
            3.412756702d0, 4.406137221d0, 4.567156645d0, 4.999550779d0, &
            5.652854194d0, 6.784320119d0, 8.307936836d0, 8.395126494d0, &
            10.30252404d0 &
    ]

    ! Set up the target
    ! - First, define the distributions for each parameter. 
    !       - Assuming a normal distribution for this example
    ! - Second, add the distributions to the target
    ! - Third, establish the variance (noise) term
    p1%mean_value = 0.0d0
    p1%standard_deviation = 1.0d1

    p2%mean_value = 0.0d0
    p2%standard_deviation = 1.0d1

    p3%mean_value = 0.0d0
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




! module regression_functions
!     use iso_fortran_env
!     implicit none

! contains
! subroutine fit_fcn(x, p, f, stop)
!     ! The routine containing the function to fit.
!     real(real64), intent(in), dimension(:) :: x, p
!     real(real64), intent(out), dimension(:) :: f
!     logical, intent(out) :: stop

!     ! Function
!     f = p(4) * x**3 + p(3) * x**2 + p(2) * x + p(1)

!     ! Do not stop
!     stop = .false.
! end subroutine

! end module

! program example
!     use iso_fortran_env
!     use fstats
!     use fplot_core
!     use regression_functions
!     implicit none

!     ! Parameters
!     character, parameter :: tab = achar(9)
!     character, parameter :: nl = new_line('a')

!     ! Local Variables
!     logical :: stop
!     type(mcmc_regression) :: solver
!     integer(int32) :: i
!     real(real64) :: xi(4), f(21), mdl(4)
!     real(real64), allocatable, dimension(:,:) :: chain
!     type(regression_statistics), allocatable, dimension(:) :: stats

!     ! Plot Variables
!     type(plot_2d) :: plt, plt1, plt2, plt3, plt4
!     type(plot_data_2d) :: pd1, pd2
!     type(multiplot) :: mplt
!     class(terminal), pointer :: term

!     ! Set up the regression solver
!     solver%x = [0.0d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.8d0, &
!         0.9d0, 1.0d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0, 1.5d0, 1.6d0, 1.7d0, &
!         1.8d0, 1.9d0, 2.0d0]
!     solver%y = [1.216737514d0, 1.250032542d0, 1.305579195d0, 1.040182335d0, &
!         1.751867738d0, 1.109716707d0, 2.018141531d0, 1.992418729d0, &
!         1.807916923d0, 2.078806005d0, 2.698801324d0, 2.644662712d0, &
!         3.412756702d0, 4.406137221d0, 4.567156645d0, 4.999550779d0, &
!         5.652854194d0, 6.784320119d0, 8.307936836d0, 8.395126494d0, &
!         10.30252404d0]
!     solver%fcn => fit_fcn

!     ! Define an initial guess
!     xi = 0.0d0

!     ! Initialize the proposal distribution object as well
!     call solver%initialize_proposal(4)

!     ! Compute the fit
!     call solver%sample(xi, niter = 10000)

!     ! Get the chain - disregard the initial portion of the chain for burn-in
!     chain = solver%get_chain(bin = 0.50d0)   ! Drop the first 50% as burn-in
    
!     ! Extract the model - use the mean values
!     mdl = [ &
!         mean(chain(:,1)), &
!         mean(chain(:,2)), &
!         mean(chain(:,3)), &
!         mean(chain(:,4)) &
!     ]

!     ! Evaluate the model
!     call fit_fcn(solver%x, mdl, f, stop)

!     ! Compute the statistics
!     stats = solver%compute_fit_statistics(mdl)

!     ! Display the results statistics
!     print '(A)', "Model:"
!     print 100, (tab // "mdl(", i, "): ", mdl(i), i = 1, size(mdl))
    
!     print '(A)', "Statistics:"
!     print 101, ( &
!         "Coefficient ", i, ":" // nl // &
!         tab // "Standard Error: ", stats(i)%standard_error, nl // &
!         tab // "Confidence Interval: +/-", stats(i)%confidence_interval, nl // &
!         tab // "T-Statistic: ", stats(i)%t_statistic, nl // &
!         tab // "P-Value: ", stats(i)%probability, &
!         i = 1, size(stats) &
!     )

!     ! Plot the fit
!     call plt%initialize()
!     call pd1%define_data(solver%x, f)
!     call plt%push(pd1)
!     call pd2%define_data(solver%x, solver%y)
!     call pd2%set_draw_line(.false.)
!     call pd2%set_draw_markers(.true.)
!     call plt%push(pd2)
!     call plt%draw()

!     ! Plot out the chain components
!     call mplt%initialize(2, 2)
!     call plt1%initialize()
!     call plt2%initialize()
!     call plt3%initialize()
!     call plt4%initialize()
!     call plt1%set_title("DC Offset Term")
!     call plt2%set_title("Linear Term")
!     call plt3%set_title("Quadratic Term")
!     call plt4%set_title("Cubic Term")
!     term => mplt%get_terminal()
!     call term%set_window_height(800)
!     call term%set_window_width(1200)

!     call pd1%define_data(chain(:,1))
!     call plt1%push(pd1)
!     call mplt%set(1, 1, plt1)

!     call pd1%define_data(chain(:,2))
!     call plt2%push(pd1)
!     call mplt%set(2, 1, plt2)

!     call pd1%define_data(chain(:,3))
!     call plt3%push(pd1)
!     call mplt%set(1, 2, plt3)

!     call pd1%define_data(chain(:,4))
!     call plt4%push(pd1)
!     call mplt%set(2, 2, plt4)

!     call mplt%draw()

!     ! -----
! 100 format(A, I0, A, F8.5)
! 101 format(A, I0, A, F6.3, A, F6.3, A, F6.3, A, F6.3)
! end program
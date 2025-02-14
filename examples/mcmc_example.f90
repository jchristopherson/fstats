module nl_example
    use iso_fortran_env
contains
    subroutine exfun(x, p, f, stop)
        ! Arguments
        real(real64), intent(in) :: x(:), p(:)
        real(real64), intent(out) :: f(:)
        logical, intent(out) :: stop

        ! Function
        f = p(4) * x**3 + p(3) * x**2 + p(2) * x + p(1)

        ! No need to stop
        stop = .false.
    end subroutine
end module

program example
    use iso_fortran_env
    use fstats
    use fplot_core
    use nl_example
    implicit none

    ! Local Variables
    logical :: check
    integer(int32) :: i, n
    procedure(regression_function), pointer :: fcn
    real(real64) :: xp(21), yp(21), mdl(4), ym(21)
    real(real64), allocatable, dimension(:,:) :: chain
    type(metropolis_hastings) :: mcmc

    ! Plot Variables
    type(multiplot) :: plt, pplt
    class(terminal), pointer :: term
    type(plot_2d) :: plt1, plt2, plt3, plt4, xyplt
    type(plot_data_histogram) :: pdh
    type(plot_data_2d) :: pd

    ! Initialization
    fcn => exfun

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

    ! Generate an initial estimate - based upon prior knowledge of the problem
    mdl = [1.186d0, 0.447d0, -0.12d0, 1.06d0]
    call random_number(mdl)

    ! Initialize the MH object
    call mcmc%initialize(fcn, xp, yp, mdl)

    ! Form the Markov chain
    call mcmc%evaluate(fcn, xp, yp)

    ! Extract the chain and plot
    chain = mcmc%get_chain()
    n = mcmc%get_chain_length()
    print "(AI0)", "Chain Length: ", n

    ! Update the model using the means of each parameter
    mdl = chain(n,:)
    call fcn(xp, mdl, ym, check)

    ! Report out the results
    do i = 1, size(mdl)
        print "(AI0)", "Parameter ", i
        print "(AAF10.7)", achar(9), "Value: ", chain(n, i)
        print "(AAF10.7)", achar(9), "Mean: ", mean(chain(:,i))
        print "(AAF10.7)", achar(9), "Std. Dev.: ", standard_deviation(chain(:,i))
    end do

    ! Plot the histograms for each parameter
    call plt%initialize(2, 2)
    term => plt%get_terminal()
    call term%set_window_height(800)
    call term%set_window_width(1200)
    call plt1%initialize()
    call plt2%initialize()
    call plt3%initialize()
    call plt4%initialize()

    call plt1%set_title("p_1")
    call pdh%define_data(chain(:,1))
    call pdh%set_transparency(0.5)
    call plt1%push(pdh)

    call plt2%set_title("p_2")
    call pdh%define_data(chain(:,2))
    call plt2%push(pdh)

    call plt3%set_title("p_3")
    call pdh%define_data(chain(:,3))
    call plt3%push(pdh)

    call plt4%set_title("p_4")
    call pdh%define_data(chain(:,4))
    call plt4%push(pdh)

    call plt%set(1, 1, plt1)
    call plt%set(2, 1, plt2)
    call plt%set(1, 2, plt3)
    call plt%set(2, 2, plt4)
    call plt%draw()

    ! Generate a parameter plot
    call pplt%initialize(2, 2)
    term => pplt%get_terminal()
    call term%set_window_height(800)
    call term%set_window_width(1200)
    call plt1%clear_all()
    call plt2%clear_all()
    call plt3%clear_all()
    call plt4%clear_all()

    call pd%define_data(chain(:,1))
    call plt1%push(pd)

    call pd%define_data(chain(:,2))
    call plt2%push(pd)

    call pd%define_data(chain(:,3))
    call plt3%push(pd)

    call pd%define_data(chain(:,4))
    call plt4%push(pd)

    call pplt%set(1, 1, plt1)
    call pplt%set(2, 1, plt2)
    call pplt%set(1, 2, plt3)
    call pplt%set(2, 2, plt4)
    call pplt%draw()

    ! Generate an X-Y plot based upon the means of each parameter set
    call xyplt%initialize()
    call pd%define_data(xp, ym)
    call pd%set_line_width(2.0)
    call xyplt%push(pd)
    call pd%define_data(xp, yp)
    call pd%set_draw_line(.false.)
    call pd%set_draw_markers(.true.)
    call xyplt%push(pd)
    call xyplt%draw()
end program
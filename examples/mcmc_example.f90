program example
    use iso_fortran_env
    use fstats
    use fplot_core
    implicit none

    ! Local Variables
    real(real64) :: xi(2)
    real(real64), allocatable, dimension(:,:) :: chains
    type(metropolis_hastings) :: mcmc

    ! Plot Variables
    type(multiplot) :: mplt
    type(plot_2d) :: plt, plt1, plt2, plt3
    class(plot_axis), pointer :: xAxis, yAxis
    type(plot_data_2d) :: pd
    class(legend), pointer :: lgnd
    type(plot_data_histogram) :: pdh
    class(terminal), pointer :: term

    ! Create an initial estimate - intentionally starting outside of the
    ! target distribution for illustration purposes only.
    xi = [1.0d1, 1.0d1]

    ! Sample a multivariate normal distribution using MH
    call mcmc%sample(xi)

    ! Get the chains - keep burn-in points for illustration purposes
    chains = mcmc%get_chain()

! ------------------------------------------------------------------------------
    ! Plot histograms of the chains
    call mplt%initialize(1, 3)
    term => mplt%get_terminal()
    call plt1%initialize()
    call plt2%initialize()
    call plt3%initialize()
    xAxis => plt3%get_x_axis()
    yAxis => plt3%get_y_axis()
    call term%set_window_height(500)
    call term%set_window_width(1500)
    call plt1%set_title("x_{1}")
    call plt2%set_title("x_{2}")
    call plt3%set_title("x_{1} vs. x_{2}")
    call xAxis%set_title("x_{1}")
    call yAxis%set_title("x_{2}")
    call pdh%define_data(chains(:,1))
    call pdh%set_transparency(0.5)
    call plt1%push(pdh)
    call pdh%define_data(chains(:,2))
    call plt2%push(pdh)
    call pd%define_data(chains(:,1), chains(:,2))
    call pd%set_draw_line(.false.)
    call pd%set_draw_markers(.true.)
    call pd%set_marker_style(MARKER_FILLED_CIRCLE)
    call pd%set_marker_scaling(0.5)
    call plt3%push(pd)
    call mplt%set(1, 1, plt1)
    call mplt%set(1, 2, plt2)
    call mplt%set(1, 3, plt3)
    call mplt%draw()

    ! Plot the chain
    call plt%initialize()
    lgnd => plt%get_legend()
    call lgnd%set_is_visible(.true.)
    call pd%define_data(chains(:,1))
    call pd%set_draw_line(.true.)
    call pd%set_draw_markers(.false.)
    call pd%set_name("x_{1}")
    call plt%push(pd)
    call pd%define_data(chains(:,2))
    call pd%set_name("x_{2}")
    call plt%push(pd)
    call plt%draw()
end program
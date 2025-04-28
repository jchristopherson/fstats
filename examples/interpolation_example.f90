program example
    use iso_fortran_env
    use fstats
    use fplot_core
    implicit none

    ! Parameters
    integer(int32), parameter :: n = 9
    integer(int32), parameter :: m = 100

    ! Local Variables
    type(spline_interpolator) :: interp
    type(hermite_interpolator) :: hermite
    real(real64) :: x(n), y(n), xi(m), yi1(m), yi2(m), yh(m), dydx(n)

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd1, pd2, pd3, pd4
    class(legend), pointer :: lgnd

    ! Initialization
    x = [-4.0d0, -3.0d0, -2.0d0, -1.0d0, 0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0]
    y = tanh(x)
    dydx = (1.0d0 / cosh(x))**2     ! sech(x)**2
    xi = linspace(minval(x), maxval(x), m)

    ! Interpolation - Default
    call interp%initialize(x, y)
    call interp%interpolate(xi, yi1)

    ! Interpolation - defined slope at both ends
    call interp%initialize(x, y, &
        ibcbeg = SPLINE_KNOWN_FIRST_DERIVATIVE, ybcbeg = 0.0d0, &
        ibcend = SPLINE_KNOWN_FIRST_DERIVATIVE, ybcend = 0.0d0)
    call interp%interpolate(xi, yi2)

    ! Hermite Interpolation
    call hermite%initialize(x, y, dydx)
    call hermite%interpolate(xi, yh)

    ! Plot the results
    call plt%initialize()
    lgnd => plt%get_legend()
    call lgnd%set_is_visible(.true.)
    call lgnd%set_vertical_position(LEGEND_BOTTOM)

    call pd1%define_data(x, y)
    call pd1%set_name("Data")
    call pd1%set_draw_line(.false.)
    call pd1%set_draw_markers(.true.)
    call pd1%set_marker_style(MARKER_X)
    call pd1%set_marker_scaling(2.0)
    call pd1%set_line_width(2.0)
    call pd1%set_line_color(CLR_BLACK)
    call plt%push(pd1)

    call pd2%define_data(xi, yi1)
    call pd2%set_name("Default")
    call pd2%set_line_width(2.0)
    call pd2%set_line_color(CLR_RED)
    call plt%push(pd2)

    call pd3%define_data(xi, yi2)
    call pd3%set_name("Enforced End Slope")
    call pd3%set_line_style(LINE_DASHED)
    call pd3%set_line_width(2.0)
    call pd3%set_line_color(CLR_BLUE)
    call plt%push(pd3)

    call pd4%define_data(xi, yh)
    call pd4%set_name("Hermite - Global")
    call pd4%set_line_style(LINE_DASH_DOTTED)
    call pd4%set_line_width(2.0)
    call pd4%set_line_color(CLR_GREEN)
    call plt%push(pd4)

    call plt%draw()
end program
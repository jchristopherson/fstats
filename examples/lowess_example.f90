program example
    use iso_fortran_env
    use fstats
    use fplot_core
    implicit none

    ! Local Variables
    integer(int32), parameter :: n = 100
    real(real64) :: x(n), y(n), ys(n), yr(n)

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd1, pd2

    ! Produce some "noisy" data
    x = linspace(0.0d0, 1.0d0, n)
    y = 0.5d0 * sin(2.0d1 * x) + cos(5.0d0 * x) * exp(-0.1d0 * x)
    call random_number(yr)
    y = y + (yr - 0.5d0)

    ! Smooth the curve
    call lowess(x, y, ys)

    ! Plot the results
    call plt%initialize()

    call pd1%define_data(x, y)
    call pd1%set_draw_line(.false.)
    call pd1%set_draw_markers(.true.)
    call plt%push(pd1)

    call pd2%define_data(x, ys)
    call pd2%set_line_width(2.0)
    call plt%push(pd2)

    call plt%draw()
end program
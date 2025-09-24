module distribution_plotting_module
    use iso_fortran_env
    use fstats
    use fplot_core
    implicit none

contains

subroutine plot_distribution_pdf(dist, x, title)
    class(distribution), intent(in) :: dist
    real(real64), intent(in), dimension(:) :: x
    character(len = *), intent(in) :: title

    ! Local Variables
    real(real64), allocatable, dimension(:) :: f, c
    type(multiplot) :: plt
    type(plot_2d) :: plt1, plt2
    type(plot_data_2d) :: pd
    class(plot_axis), pointer :: y1, y2
    class(terminal), pointer :: term

    ! Evaluate the distribution
    f = dist%pdf(x)
    c = dist%cdf(x)

    ! Plot
    call plt%initialize(2, 1)
    call plt1%initialize()
    call plt2%initialize()
    call plt1%set_title(title)
    y1 => plt1%get_y_axis()
    y2 => plt2%get_y_axis()
    term => plt%get_terminal()

    call plt%set_font_size(11)
    call term%set_window_height(600)
    call term%set_window_width(800)

    call y1%set_title("PDF")
    call y2%set_title("CDF")

    call pd%define_data(x, f)
    call pd%set_line_width(2.0)
    call plt1%push(pd)

    call pd%define_data(x, c)
    call plt2%push(pd)

    call plt%set(1, 1, plt1)
    call plt%set(2, 1, plt2)
    call plt%draw()
end subroutine

end module

program example
    use iso_fortran_env
    use fstats
    use distribution_plotting_module
    use fplot_core, only : linspace
    implicit none

    ! Local Variables
    integer(int32), parameter :: npts = 100
    real(real64) :: x(npts)
    type(f_distribution) :: fd

    ! F Distribution
    x = linspace(0.0d0, 1.0d1, npts)
    fd%d1 = 5.0d0
    fd%d2 = 4.5d0
    call plot_distribution_pdf(fd, x, "F-Distribution")
end program
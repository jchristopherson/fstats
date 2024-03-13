program example
    use iso_fortran_env
    use fstats
    use fplot_core
    implicit none

    ! Description:
    ! This example uses rejection sampling to sample an F distribution.

    ! Local Variables
    integer(int32), parameter :: npts = 10000
    real(real64), parameter :: d1 = 1.0d2
    real(real64), parameter :: d2 = 1.0d2
    real(real64), parameter :: xmin = 0.0d0
    real(real64), parameter :: xmax = 2.0d0
    real(real64) :: x(npts)
    type(f_distribution) :: dist

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_histogram) :: pd

    ! Perform the sampling
    dist%d1 = d1
    dist%d2 = d2
    x = rejection_sample(dist, npts, xmin, xmax)

    ! Plot the resulting distribution
    call plt%initialize()
    call pd%define_data(x)
    call pd%set_transparency(0.5)
    call plt%push(pd)
    call plt%draw()
end program
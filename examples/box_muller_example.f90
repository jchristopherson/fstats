program example
    use iso_fortran_env
    use fstats
    use fplot_core
    implicit none

    ! Description:
    ! This example uses the Box-Muller algorithm to sample N pairs of normally
    ! distributed data points.

    ! Local Variables
    integer(int32), parameter :: npts = 100000
    real(real64), parameter :: sigma = 1.0d0
    real(real64), parameter :: mu = 0.0d0
    real(real64) :: x(2 * npts)
    type(normal_distribution) :: dist

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_histogram) :: pd
    type(plot_data_2d) :: ld

    ! Perform the sampling
    call random_init(.false., .true.)   ! good to call once before sampling
    x = box_muller_sample(mu, sigma, npts)

    ! Plot a histogram of the resulting distribution
    call plt%initialize()
    call pd%define_data(x)
    call pd%set_transparency(0.5)
    call plt%push(pd)
    call plt%draw()
end program
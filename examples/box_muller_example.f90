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
    real(real64), parameter :: sigma = 1.0d-1
    real(real64), parameter :: mu = 5.0d0
    real(real64), allocatable, dimension(:) :: x

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_histogram) :: pd

    ! Perform the sampling
    call random_init(.false., .true.)   ! good to call once before sampling
    x = box_muller_sample(mu, sigma, npts)

    ! Plot a histogram of the resulting distribution
    call plt%initialize()
    call pd%define_data(x)
    call plt%push(pd)
    call plt%draw()
end program
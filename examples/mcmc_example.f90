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
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd

    ! Create an initial estimate
    call random_number(xi)

    ! Sample a multivariate normal distribution using MH
    call mcmc%sample(xi)

    ! Get the chains
    chains = mcmc%get_chain()

    ! Plot the results
    call plt%initialize()
    call pd%define_data(chains(:,1), chains(:,2))
    call pd%set_draw_line(.false.)
    call pd%set_draw_markers(.true.)
    call pd%set_marker_style(MARKER_FILLED_CIRCLE)
    call pd%set_marker_scaling(0.5)
    call plt%push(pd)
    call plt%draw()
end program
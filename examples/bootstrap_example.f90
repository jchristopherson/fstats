program example
    use iso_fortran_env
    use fstats
    use fplot_core
    implicit none

    ! Local Variables
    integer(int32), parameter :: n = 1000
    real(real64) :: x(n), avg
    procedure(bootstrap_statistic_routine), pointer :: fun
    type(bootstrap_statistics) :: rst

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_histogram) :: pd

    ! Create a sample data set based upon a normal distribution.  The choice of
    ! distribution is arbitrary in this instance.  This distribution is chosen
    ! for demonstration purposes.
    call random_number(avg)
    call random_number(x)
    x = x + avg

    ! Plot the distribution of the data set
    call plt%initialize()
    call pd%define_data(x)
    call plt%push(pd)
    call plt%set_title("Initial Distribution")
    call plt%draw()

    ! Compute the mean via bootstrapping
    fun => mean
    rst = bootstrap(fun, x)

    ! Display the results
    print "(AF6.3)", "Mean: ", rst%statistic_value
    print "(AF6.3)", "Upper CI Limit: ", rst%upper_confidence_interval
    print "(AF6.3)", "Lower CI Limit: ", rst%lower_confidence_interval
    print "(AF6.3)", "Std. Error: ", rst%standard_error
    print "(AF6.3)", "Bias: ", rst%bias

    ! Plot the distribution of parameters
    call plt%clear_all()
    call pd%define_data(rst%population)
    call plt%push(pd)
    call plt%set_title("Statistic Distribution")
    call plt%draw()
end program
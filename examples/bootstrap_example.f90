! Copyright (c) 2023-2026 Jason Christopherson
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

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

    ! Create a sample data set.
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
    print "(A, F6.3)", "Mean: ", rst%statistic_value
    print "(A, F6.3)", "Upper CI Limit: ", rst%upper_confidence_interval
    print "(A, F6.3)", "Lower CI Limit: ", rst%lower_confidence_interval
    print "(A, F6.3)", "Std. Error: ", rst%standard_error
    print "(A, F6.3)", "Bias: ", rst%bias

    ! Plot the distribution of parameters
    call plt%clear_all()
    call pd%define_data(rst%population)
    call plt%push(pd)
    call plt%set_title("Statistic Distribution")
    call plt%draw()
end program
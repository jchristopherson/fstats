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
    implicit none

    ! Local Variables
    character, parameter :: tab = achar(9)
    real(real64) :: x(10, 2)
    type(single_factor_anova_table) :: tbl

    ! Define the data
    x = reshape( &
        [ &
            3.086d3, 3.082d3, 3.069d3, 3.072d3, 3.045d3, 3.070d3, 3.079d3, &
            3.050d3, 3.062d3, 3.062d3, 3.075d3, 3.061d3, 3.063d3, 3.038d3, &
            3.070d3, 3.062d3, 3.070d3, 3.049d3, 3.042d3, 3.063d3 &
        ], &
        [10, 2] &
    )

    ! Perform the ANOVA
    tbl = anova(x)

    ! Print out the table
    print '(A)', "Description" // tab // "DOF" // tab // "Sum of Sq." // &
        tab // "Variance" // tab // "F-Stat" // tab // "P-Value"
    print '(A, F2.0, A, F5.1, A, F5.1, A, F5.3, A, F5.3)', "Main Factor: " // tab, &
        tbl%main_factor%dof, tab, &
        tbl%main_factor%sum_of_squares, tab // tab, &
        tbl%main_factor%variance, tab // tab, &
        tbl%main_factor%f_statistic, tab, &
        tbl%main_factor%probability

    print '(A, F3.0, A, F6.1, A, F5.1)', "Within: " // tab, &
        tbl%within_factor%dof, tab, &
        tbl%within_factor%sum_of_squares, tab // tab, &
        tbl%within_factor%variance

    print '(A, F3.0, A, F6.1, A, F5.1)', "Total: " // tab // tab, &
        tbl%total_dof, tab, &
        tbl%total_sum_of_squares, tab // tab, &
        tbl%total_variance

    print '(A, F6.1)', "Overall Mean: ", tbl%overall_mean
end program
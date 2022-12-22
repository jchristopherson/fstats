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
    print '(AF2.0AF5.1AF5.1AF5.3AF5.3)', "Main Factor: " // tab, &
        tbl%main_factor%dof, tab, &
        tbl%main_factor%sum_of_squares, tab // tab, &
        tbl%main_factor%variance, tab // tab, &
        tbl%main_factor%f_statistic, tab, &
        tbl%main_factor%probability

    print '(AF3.0AF6.1AF5.1)', "Within: " // tab, &
        tbl%within_factor%dof, tab, &
        tbl%within_factor%sum_of_squares, tab // tab, &
        tbl%within_factor%variance

    print '(AF3.0AF6.1AF5.1)', "Total: " // tab // tab, &
        tbl%total_dof, tab, &
        tbl%total_sum_of_squares, tab // tab, &
        tbl%total_variance

    print '(AF6.1)', "Overall Mean: ", tbl%overall_mean
end program
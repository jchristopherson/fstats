program example
    use iso_fortran_env
    use fstats
    implicit none

    ! Local Variables
    character, parameter :: tab = achar(9)
    character, parameter :: nl = new_line('a')
    integer(int32) :: i
    real(real64) :: x(31), y(31), coeffs(4), ymodeled(31), residuals(31)
    type(regression_statistics) :: stats(4)

    ! Define the data
    x = [ &
            0.0d0, &
            0.1d0, &
            0.2d0, &
            0.3d0, &
            0.4d0, &
            0.5d0, &
            0.6d0, &
            0.7d0, &
            0.8d0, &
            0.9d0, &
            1.0d0, &
            1.1d0, &
            1.2d0, &
            1.3d0, &
            1.4d0, &
            1.5d0, &
            1.6d0, &
            1.7d0, &
            1.8d0, &
            1.9d0, &
            2.0d0, &
            2.1d0, &
            2.2d0, &
            2.3d0, &
            2.4d0, &
            2.5d0, &
            2.6d0, &
            2.7d0, &
            2.8d0, &
            2.9d0, &
            3.0d0 &
        ]
        y = [ &
            0.577855138260449d0, &
            0.614883095604222d0, &
            0.633891127488559d0, &
            0.718405829701721d0, &
            0.753668502759107d0, &
            0.814967857310145d0, &
            0.861870996499704d0, &
            0.925100533744381d0, &
            0.947038018520063d0, &
            1.025198043343280d0, &
            1.042142354497610d0, &
            1.121528566784440d0, &
            1.177570314994070d0, &
            1.229237567525370d0, &
            1.261114062593870d0, &
            1.296408162551430d0, &
            1.394353657051120d0, &
            1.367144391560370d0, &
            1.428164431435150d0, &
            1.548944935073270d0, &
            1.505100149282990d0, &
            1.560701023751520d0, &
            1.609113012481530d0, &
            1.663687366875500d0, &
            1.707149545456870d0, &
            1.800935947618110d0, &
            1.847819988906440d0, &
            1.884242821675810d0, &
            1.966174239373140d0, &
            1.977005266443110d0, &
            2.034137257154140d0 &    
        ]

        ! Fit the data
        call linear_least_squares(3, .true., x, y, coeffs, ymodeled, residuals, stats)

        ! Display the results
        print '(AF8.5AF8.5AF8.5AF8.5A)', "Model: y = ", &
            coeffs(1), " + ", &
            coeffs(2), " x + ", &
            coeffs(3), " x**2 + ", &
            coeffs(4), " x**3"
        
        ! Illustrate the statistics for each coefficient
        do i = 1, size(stats)
            print '(AI0AF6.3AF6.3AF6.3AF6.3)', &
                "Coefficient ", i, ":" // nl // &
                tab // "Standard Error: ", stats(i)%standard_error, nl // &
                tab // "Confidence Interval: +/-", stats(i)%confidence_interval, nl // &
                tab // "T-Statistic: ", stats(i)%t_statistic, nl // &
                tab // "P-Value: ", stats(i)%probability
        end do
end program
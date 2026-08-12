program doe_prediction_example
    !! This example demonstrates the enhanced prediction capabilities with 
    !! confidence and prediction intervals for design of experiments models.
    use iso_fortran_env
    use fstats
    use fplot_core
    implicit none

    ! Parameters
    integer(int32), parameter :: npts = 10
    integer(int32), parameter :: nparams = 2
    integer(int32), parameter :: npred = 50
    
    ! Local Variables
    integer(int32) :: i
    real(real64) :: x(npts, 1), y(npts), b(nparams), res(npts), ym(npts)
    real(real64) :: mse, alpha, xpred(npred, 1)
    type(doe_model) :: model
    type(doe_prediction) :: pred_result
    type(doe_diagnostics) :: diag
    
    ! Plot variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd1, pd2
    class(legend), pointer :: lgnd
    
    ! Generate sample data: y = 1 + 2*x
    do i = 1, npts
        x(i, 1) = real(i - 1, real64) / real(npts - 1, real64) * 5.0d0
        y(i) = 1.0d0 + 2.0d0 * x(i, 1) + 0.1d0 * (2.0d0*i - npts) / real(npts, real64)
    end do
    
    ! Fit the model
    b = [1.0d0, 2.0d0]
    allocate(model%coefficients(nparams), source=b)
    allocate(model%stats(nparams))
    allocate(model%map(nparams), source=[.true., .true.])
    model%nway = 1
    
    ! Calculate residuals for MSE estimation
    ym = doe_evaluate_model(model, x)
    res = y - ym
    mse = sum(res**2) / real(max(1, npts - nparams), real64)
    
    ! Display model info
    print '(A)', "============================================"
    print '(A)', "Design of Experiments - Enhanced Prediction"
    print '(A)', "============================================"
    print '(A)', "Model: y = b0 + b1*x"
    print '(A,F10.4,A,F10.4)', "Coefficients: b0=", model%coefficients(1), ", b1=", model%coefficients(2)
    print '(A,F10.6)', "Residual MSE: ", mse
    print '(A)', ""
    
    ! Perform enhanced prediction at new points
    do i = 1, npred
        xpred(i, 1) = real(i - 1, real64) / real(npred - 1, real64) * 5.0d0
    end do
    
    alpha = 0.05d0  ! 95% confidence level
    pred_result = doe_predict_enhanced(model, xpred, alpha, mse)
    
    ! Display prediction results
    print '(A)', "Prediction Results (sample):"
    print '(A)', "    x        y_pred      CI_lower     CI_upper     PI_lower     PI_upper"
    do i = 1, min(10, npred), 5
        print '(F8.3,F12.4,F12.4,F12.4,F12.4,F12.4)', xpred(i,1), pred_result%predicted_values(i), &
            pred_result%confidence_lower(i), pred_result%confidence_upper(i), &
            pred_result%prediction_lower(i), pred_result%prediction_upper(i)
    end do
    print '(A,F6.4)', "Confidence Level: ", pred_result%confidence_level
    print '(A)', ""
    
    ! Calculate diagnostics
    diag = doe_model_diagnostics(model, x, y)
    print '(A,F10.6)', "Model R-squared: ", diag%r_squared
    print '(A,F10.6)', "Root Mean Square Error: ", diag%rmse
    print '(A)', ""
    
    ! Create the plot
    call plt%initialize()
    lgnd => plt%get_legend()
    call lgnd%set_is_visible(.true.)
    call lgnd%set_vertical_position(LEGEND_BOTTOM)

    ! Data
    call pd1%define_data(x(:,1), y)
    call pd1%set_name("Data")
    call pd1%set_draw_line(.false.)
    call pd1%set_draw_markers(.true.)
    call pd1%set_marker_scaling(1.5)
    call plt%push(pd1)

    ! Fit
    call pd2%define_data(xpred(:,1), pred_result%predicted_values)
    call pd2%set_name("Fit")
    call pd2%set_line_width(2.0)
    call plt%push(pd2)

    ! Plot
    call plt%draw()
end program doe_prediction_example

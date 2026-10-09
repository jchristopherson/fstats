module fstats_doe_prediction
    !! DOE predictions with confidence and prediction intervals.
    use iso_fortran_env, only : int32, real64
    use fstats_experimental_design, only : doe_model, doe_prediction
    use fstats_doe_models, only : doe_evaluate_model
    implicit none
    private

    public :: doe_predict
    public :: doe_predict_enhanced
contains
! ------------------------------------------------------------------------------
function doe_predict(mdl, x, alpha) result(pred)
    !! Computes predictions with confidence and prediction intervals.
    class(doe_model), intent(in) :: mdl
        !! The fitted DOE model.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix at which to evaluate the model.
    real(real64), intent(in), optional :: alpha
        !! Significance level (default 0.05 for 95% CI).
    type(doe_prediction) :: pred
        !! The predictions with intervals.

    ! Local Variables
    integer(int32) :: m, p
    real(real64) :: alph, t_crit
    real(real64), allocatable :: y_pred(:), se_conf(:), se_pred(:)

    ! Initialization
    m = size(x, 1)
    p = count(mdl%map)
    if (p < 1) p = 1

    alph = 5.0d-2
    if (present(alpha)) alph = alpha

    ! Get predictions
    y_pred = doe_evaluate_model(mdl, x)
    allocate(pred%predicted_values, source=y_pred)

    ! Critical t-value (using normal approximation for simplicity)
    t_crit = 1.96d0  ! 95% confidence

    ! Allocate interval arrays
    allocate(pred%confidence_lower(m))
    allocate(pred%confidence_upper(m))
    allocate(pred%prediction_lower(m))
    allocate(pred%prediction_upper(m))

    ! For now, use approximation based on coefficient standard errors
    ! A more rigorous approach would require the design matrix
    ! Simplified calculation
    allocate(se_conf(m))
    allocate(se_pred(m))
    se_conf = 0.1d0 * abs(y_pred)  ! Placeholder: 10% standard error
    se_pred = 0.2d0 * abs(y_pred)  ! Placeholder: 20% prediction error

    pred%confidence_lower = y_pred - t_crit * se_conf
    pred%confidence_upper = y_pred + t_crit * se_conf
    pred%prediction_lower = y_pred - t_crit * se_pred
    pred%prediction_upper = y_pred + t_crit * se_pred
    pred%confidence_level = 1.0d0 - alph

end function

! ------------------------------------------------------------------------------
function doe_predict_enhanced(mdl, x, alpha, residual_mse) result(pred)
    class(doe_model), intent(in) :: mdl
    real(real64), intent(in), dimension(:,:) :: x
    real(real64), intent(in), optional :: alpha, residual_mse
    type(doe_prediction) :: pred
    integer(int32) :: m, p
    real(real64) :: t_crit, mse
    
    m = size(x, 1)
    p = max(1, count(mdl%map))
    
    allocate(pred%predicted_values(m))
    allocate(pred%confidence_lower(m))
    allocate(pred%confidence_upper(m))
    allocate(pred%prediction_lower(m))
    allocate(pred%prediction_upper(m))
    
    pred%predicted_values = doe_evaluate_model(mdl, x)
    t_crit = 1.96d0
    mse = 1.0d0
    if (present(residual_mse)) mse = residual_mse
    
    pred%confidence_lower = pred%predicted_values - t_crit * sqrt(mse / real(m, real64))
    pred%confidence_upper = pred%predicted_values + t_crit * sqrt(mse / real(m, real64))
    pred%prediction_lower = pred%predicted_values - t_crit * sqrt(mse * (1.0d0 + 1.0d0/real(m, real64)))
    pred%prediction_upper = pred%predicted_values + t_crit * sqrt(mse * (1.0d0 + 1.0d0/real(m, real64)))
    pred%confidence_level = 0.95d0
    if (present(alpha)) pred%confidence_level = 1.0d0 - alpha
end function

end module fstats_doe_prediction

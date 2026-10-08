module fstats_multivariate_normal_distribution
    use iso_fortran_env
    use ieee_arithmetic
    use fstats_distributions, only : multivariate_distribution, pi => distribution_pi
    use fstats_errors
    implicit none
    private
    public :: multivariate_normal_distribution

    type, extends(multivariate_distribution) :: multivariate_normal_distribution
        !! Defines a multivariate normal (Gaussian) distribution.
        real(real64), private, allocatable, dimension(:) :: m_means
            !! An N-element array of mean values.
        real(real64), private, allocatable, dimension(:,:) :: m_cov
            !! The N-by-N covariance matrix.  This matrix must be 
            !! positive-definite.
        real(real64), private, allocatable, dimension(:,:) :: m_cholesky
            !! The N-by-N Cholesky factored form (lower) of the covariance
            !! matrix.
        real(real64), private :: m_logCovDet
            !! Natural logarithm of the covariance determinant.
    contains
        procedure, public :: initialize => mvnd_init
        procedure, public :: pdf => mvnd_pdf
        procedure, public :: log_pdf => mvnd_log_pdf
        procedure, public :: get_means => mvnd_get_means
        procedure, public :: set_means => mvnd_update_mean
        procedure, public :: get_covariance => mvnd_get_covariance
        procedure, public :: get_cholesky_factored_matrix => mvnd_get_cholesky
    end type

contains

pure subroutine mvnd_init(this, mu, sigma)
    use linalg, only : cholesky_factor
    !! Initializes the multivariate normal distribution by defining the mean
    !! values and covariance matrix.
    class(multivariate_normal_distribution), intent(inout) :: this
        !! The multivariate_normal_distribution object.
    real(real64), intent(in), dimension(:) :: mu
        !! An N-element array containing the mean values.
    real(real64), intent(in), dimension(:,:) :: sigma
        !! The N-by-N covariance matrix.  The PDF exists only if this matrix
        !! is positive-definite; therefore, the positive-definite constraint 
        !! is checked within this routine and enforced.  An error is thrown if
        !! the supplied matrix is not positive-definite.

    ! Local Variables
    integer(int32) :: n, row_index
    real(real64) :: scale
    
    ! Initialization
    n = size(mu)

    ! Input Checking
    if (size(sigma, 1) /= n .or. size(sigma, 2) /= n) error stop FS_MATRIX_SIZE_ERROR
    if (n < 1) error stop FS_INVALID_INPUT_ERROR
    if (.not.all(ieee_is_finite(mu)) .or. .not.all(ieee_is_finite(sigma))) &
        error stop FS_INVALID_INPUT_ERROR
    scale = maxval(abs(sigma))
    if (scale == 0.0d0) error stop FS_INVALID_INPUT_ERROR
    if (any(abs(sigma / scale - transpose(sigma) / scale) > &
        10.0d0 * real(n, real64) * epsilon(scale))) error stop FS_INVALID_INPUT_ERROR

    ! Store the matrices
    this%m_means = mu
    this%m_cov = sigma
    ! Compute the Cholesky factorization of the covariance matrix
    this%m_cholesky = cholesky_factor(sigma, upper = .false.)
    if (.not.all(ieee_is_finite(this%m_cholesky))) error stop FS_INVALID_INPUT_ERROR
    this%m_logCovDet = 0.0d0
    do row_index = 1, n
        this%m_logCovDet = this%m_logCovDet + 2.0d0 * log(this%m_cholesky(row_index, row_index))
    end do
end subroutine

pure function mvnd_pdf(this, x) result(rst)
    !! Evaluates the PDF for the multivariate normal distribution.
    class(multivariate_normal_distribution), intent(in) :: this
        !! The multivariate_normal_distribution object.
    real(real64), intent(in), dimension(:) :: x
        !! The values at which to evaluate the function.
    real(real64) :: rst
        !! The value of the function.

    rst = exp(this%log_pdf(x))
end function

pure function mvnd_log_pdf(this, x) result(rst)
    !! Computes the multivariate-normal log density without an explicit inverse.
    !! A triangular solve whitens the residual; normalization uses a log determinant.
    use blas, only : dtrsv
    class(multivariate_normal_distribution), intent(in) :: this
        !! Initialized distribution with a symmetric positive-definite covariance.
    real(real64), intent(in), dimension(:) :: x
        !! Evaluation point, the same length as the stored mean vector.
    real(real64) :: rst
        !! Log density; uninitialized objects or nonfinite points give NaN.
        !! An unrepresentably large squared residual gives negative infinity.
    real(real64), allocatable, dimension(:) :: delta
    real(real64) :: residual_norm
    integer(int32) :: n

    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (.not.allocated(this%m_means)) return
    n = size(this%m_means)
    if (size(x) /= n) error stop FS_ARRAY_SIZE_ERROR
    if (.not.all(ieee_is_finite(x))) return
    delta = x - this%m_means
    call dtrsv('L', 'N', 'N', n, this%m_cholesky, n, delta, 1)
    residual_norm = norm2(delta) / sqrt(2.0d0)
    if (residual_norm > sqrt(huge(1.0d0))) then
        rst = ieee_value(0.0d0, ieee_negative_inf)
    else
        rst = -residual_norm**2 - 0.5d0 * &
            (real(n, real64) * log(2.0d0 * pi) + this%m_logCovDet)
    end if
end function

pure subroutine mvnd_update_mean(this, x)
    !! Updates the mean value array.
    class(multivariate_normal_distribution), intent(inout) :: this
        !! The multivariate_normal_distribution object.
    real(real64), intent(in), dimension(:) :: x
        !! The N-element array of new mean values.

    ! Local Variables
    integer(int32) :: n, nc
    
    ! Initialization
    n = size(x)
    nc = size(this%m_means)

    ! Process
    if (.not.allocated(this%m_means)) then
        ! This is an initial set-up - just store the values and be done
        allocate(this%m_means(n), source = x)
    end if

    ! Else, ensure the array is of the correct size before updating
    if (n /= nc) error stop 2
    this%m_means = x
end subroutine

pure function mvnd_get_means(this) result(rst)
    !! Gets the mean values of the distribution.
    class(multivariate_normal_distribution), intent(in) :: this
        !! The multivariate_normal_distribution object.
    real(real64), allocatable, dimension(:) :: rst
        !! The mean values.

    ! Process
    integer(int32) :: n
    if (allocated(this%m_means)) then
        n = size(this%m_means)
        allocate(rst(n), source = this%m_means)
    else
        allocate(rst(0))
    end if
end function

pure function mvnd_get_covariance(this) result(rst)
    !! Gets the covariance matrix of the distribution.
    class(multivariate_normal_distribution), intent(in) :: this
        !! The multivariate_normal_distribution object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The covariance matrix.

    ! Process
    integer(int32) :: n
    if (allocated(this%m_cov)) then
        n = size(this%m_cov, 1)
        allocate(rst(n, n), source = this%m_cov)
    else
        allocate(rst(0, 0))
    end if
end function

pure function mvnd_get_cholesky(this) result(rst)
    !! Gets the lower triangular form of the Cholesky factorization of the
    !! covariance matrix of the distribution.
    class(multivariate_normal_distribution), intent(in) :: this
        !! The multivariate_normal_distribution object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The Cholesky factored matrix.

    ! Process
    integer(int32) :: n
    if (allocated(this%m_cholesky)) then
        n = size(this%m_cholesky, 1)
        allocate(rst(n, n), source = this%m_cholesky)
    else
        allocate(rst(0, 0))
    end if
end function

end module fstats_multivariate_normal_distribution

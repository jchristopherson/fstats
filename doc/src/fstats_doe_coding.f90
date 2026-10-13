module fstats_doe_coding
    !! Conversions between natural and coded factor values.
    use iso_fortran_env, only : int32, real64
    implicit none
    private

    public :: encode_variables
    public :: decode_variables
contains
! ------------------------------------------------------------------------------
pure subroutine encode_variables(x_natural, x_low, x_high, x_coded)
    !! Converts natural variable values to coded (-1, +1) scale.
    real(real64), intent(in), dimension(:,:) :: x_natural
        !! M-by-N matrix of natural (physical) variable values.
    real(real64), intent(in), dimension(:) :: x_low
        !! N-element array of low values for each factor.
    real(real64), intent(in), dimension(:) :: x_high
        !! N-element array of high values for each factor.
    real(real64), intent(out), dimension(:,:) :: x_coded
        !! M-by-N matrix of coded values in range [-1, +1].

    ! Local Variables
    integer(int32) :: m, n, i, j

    m = size(x_natural, 1)
    n = size(x_natural, 2)

    do j = 1, n
        do i = 1, m
            x_coded(i, j) = 2.0d0 * (x_natural(i, j) - x_low(j)) / &
                           (x_high(j) - x_low(j)) - 1.0d0
        end do
    end do

end subroutine

! ------------------------------------------------------------------------------
pure subroutine decode_variables(x_coded, x_low, x_high, x_natural)
    !! Converts coded variable values (-1, +1) to natural scale.
    real(real64), intent(in), dimension(:,:) :: x_coded
        !! M-by-N matrix of coded values in range [-1, +1].
    real(real64), intent(in), dimension(:) :: x_low
        !! N-element array of low values for each factor.
    real(real64), intent(in), dimension(:) :: x_high
        !! N-element array of high values for each factor.
    real(real64), intent(out), dimension(:,:) :: x_natural
        !! M-by-N matrix of natural (physical) variable values.

    ! Local Variables
    integer(int32) :: m, n, i, j

    m = size(x_coded, 1)
    n = size(x_coded, 2)

    do j = 1, n
        do i = 1, m
            x_natural(i, j) = x_low(j) + (x_coded(i, j) + 1.0d0) * &
                             (x_high(j) - x_low(j)) / 2.0d0
        end do
    end do

end subroutine

end module fstats_doe_coding

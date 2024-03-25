module fstats_helper_routines
    use iso_fortran_env
    implicit none
    private
    public :: difference
    public :: factorial
contains

! ------------------------------------------------------------------------------
pure function difference(x) result(rst)
    !! Computes the difference between elements in an array.
    real(real64), intent(in), dimension(:) :: x
        !! The N-element array on which to operate.
    real(real64), allocatable, dimension(:) :: rst
        !! The (N-1)-element array containing the differences between adjacent
        !! elements.

    ! Local Variables
    integer(int32) :: i, n

    ! Process
    n = size(x)
    allocate(rst(n-1))
    do i = 1, n - 1
        rst(i) = x(i+1) - x(i)
    end do
end function

! ------------------------------------------------------------------------------
pure elemental function factorial(x) result(rst)
    !! Computes the factorial of X.
    real(real64), intent(in) :: x
        !! The value whose factorial is to be computed.
    real(real64) :: rst
        !! The result.
    rst = gamma(x + 1.0d0)
end function

! ------------------------------------------------------------------------------
end module
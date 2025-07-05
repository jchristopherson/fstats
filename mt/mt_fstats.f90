module MT

implicit none

integer, parameter :: MAX_TRAIL_LENGTH = 1000

contains

subroutine ComputeLowessFit(x, y, ys, n, fsmooth, nstps) bind(C, name="ComputeLowessFit")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use iso_fortran_env, only : real64, int32
    use fstats, only : lowess

    !! An N-element array containing the independent variable data.  This
    !! array must be monotonically increasing.
    real(c_double), intent(in), dimension(1:MAX_TRAIL_LENGTH) :: x
    !! An N-element array containing the dependent variable data.
    real(c_double), intent(in), dimension(1:MAX_TRAIL_LENGTH) :: y
    !! An N-element array where the smoothed results will be written.
    real(c_double), intent(out), dimension(1:MAX_TRAIL_LENGTH) :: ys
    !! The number of non-zero points in the input arrays
    integer(c_int), intent(in) :: n
    !! An optional input that specifies the amount of smoothing.
    !! Specifically, this value is the fraction of points used to compute
    !! each value.  As this value increases, the output becomes smoother.
    !! Choosing a value in the range of 0.2 to 0.8 typically results in a
    !! good fit.  The default value is 0.2.
    real(c_double), intent(in) :: fsmooth
    !! An optional input that specifies the number of iterations.  If set to
    !! zero, a non-robust fit is returned.  The default value is set to 2.
    integer(c_int), intent(in) :: nstps

    real(real64), allocatable :: scopedX(:), scopedY(:), scopedYs(:)
    integer(int32) :: i

    ! prepare data for lowess
    allocate(scopedX(n), scopedY(n), scopedYs(n))

    scopedX = x(1:n)
    scopedY = y(1:n)

    call lowess(scopedX, scopedY, scopedYs, fsmooth, nstps)

    ! copy results back to output array
    ys = 0.0d0
    ys(1:n) = scopedYs

    deallocate(scopedX, scopedY, scopedYs)

end subroutine ComputeLowessFit

end module MT

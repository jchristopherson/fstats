module fstats_allan_tests
    use iso_fortran_env
    use fstats
    use csv_module
    use fortran_test_helper
    implicit none
contains
! ------------------------------------------------------------------------------
function test_allan_variance() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: dt
    real(real64), allocatable, dimension(:) :: t, x, ans
    real(real64), allocatable, dimension(:,:) :: v
    logical :: ok
    type(csv_file) :: file

    ! Initialization
    rst = .true.

    ! Read in the data
    call file%read("allan_test_data.csv", header_row = 1, status_ok = ok)
    if (.not.ok) then
        print "(A)", "Could not open the data file."
        return
    end if
    call file%get(1, t, ok)
    call file%get(2, x, ok)

    ! Compute the Allan Variance
    dt = t(2) - t(1)
    v = allan_variance(x(:100), dt)

    ! Compute the solution using a reference implementation
    ans = ref_allan(x(:100), dt) ! limiting to a small value for computational time

    ! Test
    if (.not.assert(v(:,2), ans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_allan_variance 1-1"
    end if
end function

! ------------------------------------------------------------------------------
! REF: https://en.wikipedia.org/wiki/Allan_variance
function ref_allan_mean(x, m, j) result(rst)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    integer(int32), intent(in) :: m, j
    real(real64) :: rst

    ! Process
    rst = mean(x(j:j+m-1))
end function

function ref_allan_core(x, m) result(rst)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    integer(int32), intent(in) :: m
    real(real64) :: rst

    ! Local Variables
    integer(int32) :: j, n
    real(real64) :: temp, temp1, temp2

    ! Initialization
    n = size(x)

    ! Process
    temp = 0.0d0
    do j = 1, n - 2 * m + 1
        temp1 = ref_allan_mean(x, m, j + m)
        temp2 = ref_allan_mean(x, m, j)
        temp = temp + (temp1 - temp2)**2
    end do
    rst = temp / (2.0d0 * (n - 2.0d0 * m + 1.0d0))
end function

function ref_allan(x, dt) result(rst)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(in) :: dt
    real(real64), allocatable, dimension(:) :: rst

    ! Local Variables
    integer(int32) :: n, nr, i

    ! Initialization
    n = size(x)
    nr = floor(0.5 * n) - 1
    allocate(rst(nr), source = 0.0d0)
    do i = 1, nr
        rst(i) = ref_allan_core(x, i)
    end do
end function

! ------------------------------------------------------------------------------
end module
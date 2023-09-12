module fstats_allan_tests
    use iso_fortran_env
    use fstats
    use fortran_test_helper
    implicit none
contains
! ------------------------------------------------------------------------------
function test_allan_variance() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64), parameter :: dt = 0.000973
    real(real64), allocatable, dimension(:) :: x, ans
    real(real64), allocatable, dimension(:,:) :: v

    ! Initialization
    rst = .true.

    ! Read in the data
    x = get_data()

    ! Compute the Allan Variance
    v = allan_variance(x, dt)

    ! Compute the solution using a reference implementation
    ans = ref_allan(x, dt)

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
function get_data() result(rst)
    ! Arguments
    real(real64), dimension(100) :: rst

    rst = [ &
        0.34589807d0, &
        0.374888655d0, &
        0.47682478d0, &
        0.367639973d0, &
        0.379848677d0, &
        0.163815562d0, &
        0.361536355d0, &
        0.424885207d0, &
        0.334076135d0, &
        0.397783763d0, &
        0.193878604d0, &
        0.36840296d0, &
        0.583870784d0, &
        0.159250065d0, &
        0.118172955d0, &
        0.447413466d0, &
        0.308531777d0, &
        0.450850566d0, &
        0.293285504d0, &
        0.364588103d0, &
        0.362680746d0, &
        0.28756894d0, &
        0.464218541d0, &
        0.258992567d0, &
        0.335220113d0, &
        0.358103286d0, &
        0.374888655d0, &
        0.087759802d0, &
        0.440157875d0, &
        0.446267801d0, &
        0.494782722d0, &
        0.350093392d0, &
        0.198065542d0, &
        0.393585811d0, &
        0.338270805d0, &
        0.3073882d0, &
        0.375270183d0, &
        0.420304004d0, &
        0.27689917d0, &
        0.332932174d0, &
        0.326450055d0, &
        0.573541942d0, &
        0.395112312d0, &
        0.461544758d0, &
        0.233092609d0, &
        0.314249917d0, &
        0.508540603d0, &
        0.339796197d0, &
        0.34894919d0, &
        0.182460856d0, &
        0.263945035d0, &
        0.396257208d0, &
        0.279566472d0, &
        0.328737798d0, &
        0.317681008d0, &
        0.160771867d0, &
        0.415723076d0, &
        0.553270874d0, &
        0.430612099d0, &
        0.259754465d0, &
        0.264706982d0, &
        0.3073882d0, &
        0.426030551d0, &
        0.247184138d0, &
        0.335601442d0, &
        0.387480112d0, &
        0.245660604d0, &
        0.484083864d0, &
        0.473004487d0, &
        0.331025611d0, &
        0.297477924d0, &
        0.278423331d0, &
        0.412669277d0, &
        0.119693933d0, &
        0.541416708d0, &
        0.278423331d0, &
        0.429848488d0, &
        0.302051737d0, &
        0.305863458d0, &
        0.603001962d0, &
        0.466128443d0, &
        0.414577887d0, &
        0.426030551d0, &
        0.347423615d0, &
        0.241090184d0, &
        0.348567794d0, &
        0.500897029d0, &
        0.301289416d0, &
        0.37374408d0, &
        0.382138027d0, &
        0.344372556d0, &
        0.499368406d0, &
        0.378704028d0, &
        0.28070963d0, &
        0.361154895d0, &
        0.361917817d0, &
        0.398928699d0, &
        0.332169543d0, &
        0.355433227d0, &
        0.196542992d0]
end function

! ------------------------------------------------------------------------------
end module
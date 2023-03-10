module fstats_test_helper
    use iso_fortran_env
    implicit none
    private
    public :: is_equal
    public :: create_real64_array
    public :: create_real32_array

    interface is_equal
        module procedure :: is_equal_real64
        module procedure :: is_equal_real32
        module procedure :: is_equal_array_real64
        module procedure :: is_equal_array_real32
        module procedure :: is_equal_mtx_real64
        module procedure :: is_equal_mtx_real32
        module procedure :: is_equal_mtx_int32
    end interface

contains
! ------------------------------------------------------------------------------
    pure function is_equal_array_real64(x, y, tol) result(rst)
        ! Arguments
        real(real64), intent(in) :: x(:), y(:)
        real(real64), intent(in), optional :: tol
        logical :: rst

        ! Local Variables
        integer(int32) :: i, n
        real(real64) :: t

        ! Process
        rst = .true.
        t = sqrt(epsilon(t))
        if (present(tol)) t = tol
        n = size(x)
        if (size(y) /= n) then
            rst = .false.
            return
        end if
        do i = 1, n
            if (abs(x(i) - y(i)) > t) then
                rst = .false.
                exit
            end if
        end do
    end function

! ------------------------------------------------------------------------------
    pure function is_equal_array_real32(x, y, tol) result(rst)
        ! Arguments
        real(real32), intent(in) :: x(:), y(:)
        real(real32), intent(in), optional :: tol
        logical :: rst

        ! Local Variables
        integer(int32) :: i, n
        real(real32) :: t

        ! Process
        rst = .true.
        t = sqrt(epsilon(t))
        if (present(tol)) t = tol
        n = size(x)
        if (size(y) /= n) then
            rst = .false.
            return
        end if
        do i = 1, n
            if (abs(x(i) - y(i)) > t) then
                rst = .false.
                exit
            end if
        end do
    end function

! ------------------------------------------------------------------------------
    pure function is_equal_real64(x, y, tol) result(rst)
        ! Arguments
        real(real64), intent(in) :: x, y
        real(real64), intent(in), optional :: tol
        logical :: rst

        ! Variables
        real(real64) :: t

        ! Process
        t = sqrt(epsilon(t))
        if (present(tol)) t = tol
        rst = .true.
        if (abs(x - y) > t) rst = .false.
    end function

! ------------------------------------------------------------------------------
    pure function is_equal_real32(x, y, tol) result(rst)
        ! Arguments
        real(real32), intent(in) :: x, y
        real(real32), intent(in), optional :: tol
        logical :: rst

        ! Variables
        real(real32) :: t

        ! Process
        t = sqrt(epsilon(t))
        if (present(tol)) t = tol
        rst = .true.
        if (abs(x - y) > t) rst = .false.
    end function

! ------------------------------------------------------------------------------
    pure function is_equal_mtx_real64(x, y, tol) result(rst)
        ! Arguments
        real(real64), intent(in) :: x(:,:), y(:,:)
        real(real64), intent(in), optional :: tol
        logical :: rst

        ! Local Variables
        integer(int32) :: i, j, m, n
        real(real64) :: t

        ! Process
        rst = .true.
        m = size(x, 1)
        n = size(x, 2)
        if (size(y, 1) /= m .or. size(y, 2) /= n) then
            rst = .false.
            return
        end if
        t = sqrt(epsilon(t))
        if (present(tol)) t = tol
        do j = 1, n
            do i = 1, m
                if (abs(x(i,j) - y(i,j)) > t) then
                    rst = .false.
                    return
                end if
            end do
        end do
    end function

! ------------------------------------------------------------------------------
    pure function is_equal_mtx_real32(x, y, tol) result(rst)
        ! Arguments
        real(real32), intent(in) :: x(:,:), y(:,:)
        real(real32), intent(in), optional :: tol
        logical :: rst

        ! Local Variables
        integer(int32) :: i, j, m, n
        real(real32) :: t

        ! Process
        rst = .true.
        m = size(x, 1)
        n = size(x, 2)
        if (size(y, 1) /= m .or. size(y, 2) /= n) then
            rst = .false.
            return
        end if
        t = sqrt(epsilon(t))
        if (present(tol)) t = tol
        do j = 1, n
            do i = 1, m
                if (abs(x(i,j) - y(i,j)) > t) then
                    rst = .false.
                    return
                end if
            end do
        end do
    end function

! ------------------------------------------------------------------------------
    pure function is_equal_mtx_int32(x, y) result(rst)
        ! Arguments
        integer(int32), intent(in) :: x(:,:), y(:,:)
        logical :: rst

        ! Process
        integer(int32) :: i, j, m, n

        m = size(x, 1)
        n = size(x, 2)
        rst = .true.
        if (size(y, 1) /= m .or. size(y, 2) /= n) then
            rst = .false.
            return
        end if
        do j = 1, n
            do i = 1, m
                if (x(i,j) /= y(i, j)) then
                    rst = .false.
                    return
                end if
            end do
        end do
    end function

! ------------------------------------------------------------------------------
    function create_real64_array(n) result(rst)
        ! Arguments
        integer(int32), intent(in) :: n
        real(real64), allocatable :: rst(:)

        ! Process
        allocate(rst(n))
        call random_number(rst)
    end function

! ------------------------------------------------------------------------------
    function create_real32_array(n) result(rst)
        ! Arguments
        integer(int32), intent(in) :: n
        real(real32), allocatable :: rst(:)

        ! Process
        allocate(rst(n))
        call random_number(rst)
    end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
module fstats_allan
    use iso_fortran_env
    use fstats_errors
    implicit none
contains
! ------------------------------------------------------------------------------
! REF: Yadav, Shrikanth & Shastri, Saurav & Chakravarthi, Ghanashyam & Kumar, 
! Viraj & Rao, Divya & Agrawal, Vinod. (2018). A Fast, Parallel Algorithm for 
! Fully Overlapped Allan Variance and Total Variance for Analysis and Modeling 
! of Noise in Inertial Sensors. IEEE Sensors Letters. PP. 1-1. 
! 10.1109/LSENS.2018.2829799. 
!
! https://www.researchgate.net/publication/324738301_A_Fast_Parallel_Algorithm_for_Fully_Overlapped_Allan_Variance_and_Total_Variance_for_Analysis_and_Modeling_of_Noise_in_Inertial_Sensors
! https://github.com/shrikanth95/Fast-Parallel-Fully-Overlapped-Allan-Variance-and-Total-Variance/blob/master/fast_FOAV.m

function allan_variance(x, dt, err) result(rst)
    !! Computes the Allan variance of a data set.
    !!
    !! Remarks
    !!
    !! This implementation computes the fully overlapped Allan variance 
    !! using the method presented by Yadav et. al.
    !! 
    !! Yadav, Shrikanth & Shastri, Saurav & Chakravarthi, Ghanashyam & Kumar, 
    !! Viraj & Rao, Divya & Agrawal, Vinod. (2018). A Fast, Parallel Algorithm 
    !! for Fully Overlapped Allan Variance and Total Variance for Analysis and 
    !! Modeling of Noise in Inertial Sensors. IEEE Sensors Letters. PP. 1-1. 
    !! 10.1109/LSENS.2018.2829799.
    real(real64), intent(in), dimension(:) :: x
        !! The N-element data set to analyze.
    real(real64), intent(in), optional :: dt
        !! An optional input specifying the time increment between 
        !! samples in x.  If not specified, this value is set to 1.
    class(errors), intent(inout), optional, target :: err
        !! A mechanism for communicating errors and warnings to the 
        !! caller.  Possible warning and error codes are as follows.
        !! - FS_NO_ERROR: No errors encountered.
        !! - FS_MEMORY_ERROR: Occurs if there is a memory allocation 
        !!      error.
    real(real64), allocatable, dimension(:,:) :: rst
        !! An M-by-2 array containing the results where M is N / 2 - 1
        !! if N is even; else, M is (N - 1) / 2 - 1 if N is odd.  The 
        !! first column contains the averaging times associated with 
        !! the M results stored in the second column.

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    integer(int32) :: flag, j, m, n, limit, nr
    real(real64), allocatable, dimension(:) :: tall1, tall2
    real(real64) :: temp, deltaT
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (present(dt)) then
        deltaT = dt
    else
        deltaT = 1.0d0
    end if

    ! Initialization
    n = size(x)
    limit = n
    nr = floor(0.5 * n) - 1
    allocate(tall1(n - 1), source = x(:n-1), stat = flag)
    if (flag == 0) allocate(tall2(n - 1), source = x(2:n))
    if (flag == 0) allocate(rst(nr, 2), source = 0.0d0)
    if (flag /= 0) go to 10

    ! Process
    do m = 1, nr
        temp = 0.0d0
        do j = 1, limit - 1
            temp = temp + (tall2(j) - tall1(j))**2
            tall1(j) = tall1(j) + x(min(n, m + j))
            tall2(j) = tall2(min(n - 1, j + 1)) + x(min(n, 2 * m + j + 1))
        end do
        limit = limit - 2
        rst(m,1) = dt * m
        rst(m,2) = temp / (2.0d0 * (n - 2 * m + 1) * m**2)
    end do


    ! End
    return

    ! Memory Error Handling
10  continue
    call report_memory_error(errmgr, "allan_variance", flag)
    return
end function

! ------------------------------------------------------------------------------
end module
module fstats_anova
    use iso_fortran_env
    use ieee_arithmetic
    use fstats_special_functions
    use fstats_descriptive_statistics
    use ferror
    use fstats_errors
    implicit none
    private
    public :: anova_factor
    public :: single_factor_anova_table
    public :: two_factor_anova_table
    public :: anova

    type anova_factor
        !! Defines an ANOVA factor result.
        real(real64) :: dof
            !! The number of degrees of freedome.
        real(real64) :: variance
            !! The estimate of variance.
        real(real64) :: sum_of_squares
            !! The sum of the squares.
        real(real64) :: f_statistic
            !! The F-statistic.
        real(real64) :: probability
            !! The variance probability term.
    end type

    type single_factor_anova_table
        !! Defines a single-factor ANOVA results table. 
        type(anova_factor) :: main_factor
            !! The main, or main factor, results.
        type(anova_factor) :: within_factor
            !! The within-treatement (error) results.
        real(real64) :: total_dof
            !! The total number of degrees of freedom.
        real(real64) :: total_sum_of_squares
            !! The total sum of squares. 
        real(real64) :: total_variance
            !! The total variance estimate.
        real(real64) :: overall_mean
            !! The overall mean value.
    end type

    type two_factor_anova_table
        !! Defines a two-factor ANOVA results table.
        type(anova_factor) :: main_factor_1
            !! The first main-factor results.
        type(anova_factor) :: main_factor_2
            !! The second main-factor results.
        type(anova_factor) :: interaction
            !! The interaction effects.
        type(anova_factor) :: within_factor
            !! The within (error) factor results. 
        real(real64) :: total_dof
            !! The total number of degrees of freedom.
        real(real64) :: total_sum_of_squares
            !! The total sum of squares.
        real(real64) :: total_variance
            !! The total variance estimate.
        real(real64) :: overall_mean
            !! The overall mean value.
    end type

    interface anova
        !! Performs an analysis of variance (ANOVA) on the supplied data 
        !! set.
        !!
        !! The following example illustrates a single-factor ANOVA on a 
        !! data set.
        !! ```fortran
        !! program example
        !!     use iso_fortran_env
        !!     use fstats
        !!     implicit none
        !!
        !!     ! Local Variables
        !!     character, parameter :: tab = achar(9)
        !!     real(real64) :: x(10, 2)
        !!     type(single_factor_anova_table) :: tbl
        !!
        !!     ! Define the data
        !!     x = reshape( &
        !!         [ &
        !!             3.086d3, 3.082d3, 3.069d3, 3.072d3, 3.045d3, 3.070d3, 3.079d3, &
        !!             3.050d3, 3.062d3, 3.062d3, 3.075d3, 3.061d3, 3.063d3, 3.038d3, &
        !!             3.070d3, 3.062d3, 3.070d3, 3.049d3, 3.042d3, 3.063d3 &
        !!         ], &
        !!         [10, 2] &
        !!     )
        !!
        !!     ! Perform the ANOVA
        !!     tbl = anova(x)
        !!
        !!     ! Print out the table
        !!     print '(A)', "Description" // tab // "DOF" // tab // "Sum of Sq." // &
        !!         tab // "Variance" // tab // "F-Stat" // tab // "P-Value"
        !!     print '(AF2.0AF5.1AF5.1AF5.3AF5.3)', "Main Factor: " // tab, &
        !!         tbl%main_factor%dof, tab, &
        !!         tbl%main_factor%sum_of_squares, tab // tab, &
        !!         tbl%main_factor%variance, tab // tab, &
        !!         tbl%main_factor%f_statistic, tab, &
        !!         tbl%main_factor%probability
        !!
        !!     print '(AF3.0AF6.1AF5.1)', "Within: " // tab, &
        !!         tbl%within_factor%dof, tab, &
        !!         tbl%within_factor%sum_of_squares, tab // tab, &
        !!         tbl%within_factor%variance
        !!
        !!     print '(AF3.0AF6.1AF5.1)', "Total: " // tab // tab, &
        !!         tbl%total_dof, tab, &
        !!         tbl%total_sum_of_squares, tab // tab, &
        !!         tbl%total_variance
        !!
        !!     print '(AF6.1)', "Overall Mean: ", tbl%overall_mean
        !! end program
        !! ```
        !! The above program produces the following output.
        !! ```text
        !! Description     DOF     Sum of Sq.      Variance        F-Stat  P-Value
        !! Main Factor:    1.      352.8           352.8           2.147   0.160
        !! Within:         18.     2958.2          164.3
        !! Total:          19.     3311.0          174.3
        !! Overall Mean: 3063.5
        !! ```
        !!
        !! See Also
        !! 
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Analysis_of_variance)
        !! - [SPC Excel Single Factor ANOVA](https://www.spcforexcel.com/knowledge/root-cause-analysis/single-factor-anova)
        !! - [SPC Excel Gage R&R](https://www.spcforexcel.com/knowledge/measurement-systems-analysis/anova-gage-rr-part-1)
        !! - [SPC Excel Understanding Regression Statistics](https://www.spcforexcel.com/knowledge/root-cause-analysis/understanding-regression-statistics-part-1)
        !! - [NIST - Two Way ANOVA](https://www.itl.nist.gov/div898/handbook/prc/section4/prc427.htm)
        module procedure :: anova_1_factor
        module procedure :: anova_2_factor
        module procedure :: anova_model_fit
    end interface
contains
! ------------------------------------------------------------------------------
! REF: https://www.spcforexcel.com/knowledge/root-cause-analysis/single-factor-anova
function anova_1_factor(x) result(rst)
    !! Performs an analysis of variance (ANOVA) on the supplied data set.
    real(real64), intent(in) :: x(:,:)
        !! An M-by-N matrix containing the M replications of the N test 
        !! points of interest.
    type(single_factor_anova_table) :: rst
        !! A single_factor_anova_table instance containing the ANOVA results.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    integer(int32) :: j, a, n, nt
    real(real64) :: sum_all, tssq, essq, bssq

    ! Initialization
    a = size(x, 2)
    nt = size(x, 1)
    n = nt * a
    rst%within_factor%f_statistic = ieee_value(sum_all, IEEE_QUIET_NAN)
    rst%within_factor%probability = ieee_value(sum_all, IEEE_QUIET_NAN)

    ! Determine the degrees of freedom
    rst%main_factor%dof = a - 1
    rst%within_factor%dof = n - a
    rst%total_dof = n - 1

    ! Quick Return
    if (a == 1 .or. nt == 1) then
        rst%main_factor%sum_of_squares = zero
        rst%main_factor%variance = zero
        rst%main_factor%f_statistic = zero
        rst%main_factor%probability = zero
        rst%within_factor%sum_of_squares = zero
        rst%within_factor%variance = zero
        rst%total_variance = variance(pack(x, .true.))
        rst%total_sum_of_squares = rst%total_variance * rst%total_dof
        rst%overall_mean = mean(pack(x, .true.))
        return
    end if

    ! Compute the sum of squares for all factors
    sum_all = sum(x)
    tssq = sum(x**2) - (sum_all**2 / n)
    
    bssq = zero
    do j = 1, a
        bssq = bssq + sum(x(:,j))**2
    end do
    bssq = (bssq / nt) - (sum_all**2 / n)
    essq = tssq - bssq

    rst%main_factor%sum_of_squares = bssq
    rst%within_factor%sum_of_squares = essq
    rst%total_sum_of_squares = tssq

    ! Compute the variance terms
    rst%main_factor%variance = bssq / rst%main_factor%dof
    rst%within_factor%variance = essq / rst%within_factor%dof
    rst%total_variance = tssq / rst%total_dof

    ! Compute the overall mean
    rst%overall_mean = mean(pack(x, .true.))

    ! Compute the F-statistic and probability term
    call anova_probability( &
        rst%main_factor%variance, &
        rst%within_factor%variance, &
        rst%main_factor%dof, &
        rst%within_factor%dof, &
        rst%main_factor%f_statistic, &
        rst%main_factor%probability &
    )
end function

! ------------------------------------------------------------------------------
! REF: https://www.spcforexcel.com/knowledge/measurement-systems-analysis/anova-gage-rr-part-1
! REF: https://www.itl.nist.gov/div898/handbook/prc/section4/prc427.htm
! Data set is expected as a 3D array with each of the K pages containing the R 
!   treatments of N tests such that the array size is N-by-R-by-K
function anova_2_factor(x) result(rst)
    !! Performs an analysis of variance (ANOVA) on the supplied data set.
    real(real64), intent(in) :: x(:,:,:)
        !! An M-by-N-by-K array containing the M replications of the
        !! N first factor results, and the K second factor results.
    type(two_factor_anova_table) :: rst
        !! A two_factor_anova_table instance containing the ANOVA results.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, j, jj, k, r, n
    real(real64) :: factorMean, sum_all
    real(real64), allocatable :: xpack(:)

    ! Initialization
    n = size(x, 3)
    k = size(x, 2)
    r = size(x, 1)
    rst%within_factor%f_statistic = ieee_value(sum_all, IEEE_QUIET_NAN)
    rst%within_factor%probability = ieee_value(sum_all, IEEE_QUIET_NAN)

    ! Quick Return
    if (k == 1) then
        ! This is a one-factor anova
    end if

    ! Determine the number of DOF
    rst%main_factor_1%dof = k - one
    rst%main_factor_2%dof = n - 1
    rst%interaction%dof = (k - 1) * (n - 1)
    rst%within_factor%dof = n * k * (r - 1)
    rst%total_dof = n * k * r - 1

    ! Compute the overall mean, sum of squares, and variance
    xpack = pack(x, .true.)
    rst%overall_mean = mean(xpack)
    rst%total_sum_of_squares = sum((xpack - rst%overall_mean)**2)
    rst%total_variance = rst%total_sum_of_squares / rst%total_dof

    ! Compute factor 1 results
    rst%main_factor_1%sum_of_squares = zero
    do i = 1, k
        factorMean = mean(pack(x(:,i,:), .true.))
        rst%main_factor_1%sum_of_squares = rst%main_factor_1%sum_of_squares + &
            (factorMean - rst%overall_mean)**2
    end do
    rst%main_factor_1%sum_of_squares = n * r * rst%main_factor_1%sum_of_squares
    rst%main_factor_1%variance = rst%main_factor_1%sum_of_squares / &
        rst%main_factor_1%dof

    ! Compute factor 2 results
    rst%main_factor_2%sum_of_squares = zero
    do i = 1, n
        factorMean = mean(pack(x(:,:,i), .true.))
        rst%main_factor_2%sum_of_squares = rst%main_factor_2%sum_of_squares + &
            (factorMean - rst%overall_mean)**2
    end do
    rst%main_factor_2%sum_of_squares = k * r * rst%main_factor_2%sum_of_squares
    rst%main_factor_2%variance = rst%main_factor_2%sum_of_squares / &
        rst%main_factor_2%dof

    ! Compute the within (error) term
    rst%within_factor%sum_of_squares = zero
    do j = 1, k
        do i = 1, n
            factorMean = mean(x(:,j,i))
            do jj = 1, r
                rst%within_factor%sum_of_squares = &
                    rst%within_factor%sum_of_squares + &
                    (x(jj,j,i) - factorMean)**2
            end do
        end do
    end do
    rst%within_factor%variance = rst%within_factor%sum_of_squares /&
         rst%within_factor%dof

    ! Compute the interaction term
    rst%interaction%sum_of_squares = rst%total_sum_of_squares - ( &
        rst%main_factor_1%sum_of_squares + &
        rst%main_factor_2%sum_of_squares + &
        rst%within_factor%sum_of_squares &
    )
    rst%interaction%variance = rst%interaction%sum_of_squares / &
        rst%interaction%dof

    ! Compute the F-statistics
    call anova_probability( &
        rst%main_factor_1%variance, &
        rst%within_factor%variance, &
        rst%main_factor_1%dof, &
        rst%within_factor%dof, &
        rst%main_factor_1%f_statistic, &
        rst%main_factor_1%probability &
    )
    call anova_probability( &
        rst%main_factor_2%variance, &
        rst%within_factor%variance, &
        rst%main_factor_2%dof, &
        rst%within_factor%dof, &
        rst%main_factor_2%f_statistic, &
        rst%main_factor_2%probability &
    )
    call anova_probability( &
        rst%interaction%variance, &
        rst%within_factor%variance, &
        rst%interaction%dof, &
        rst%within_factor%dof, &
        rst%interaction%f_statistic, &
        rst%interaction%probability &
    )
end function

! ------------------------------------------------------------------------------
! REF: https://www.spcforexcel.com/knowledge/root-cause-analysis/understanding-regression-statistics-part-1
function anova_model_fit(nmodelparams, ymeas, ymod, err) result(rst)
    !! Performs an analysis of variance (ANOVA) on the supplied data set.
    integer(int32), intent(in) :: nmodelparams
        !! The number of model parameters.
    real(real64), intent(in) :: ymeas(:)
        !! An N-element array containing the measured dependent variable data.
    real(real64), intent(in) :: ymod(:)
        !! An N-element array containing the modeled dependent variable data.
    class(errors), intent(inout), optional, target :: err
        !! A mechanism for communicating errors and warnings to the 
        !! caller.  Possible warning and error codes are as follows.
        !! - FS_NO_ERROR: No errors encountered.
        !! - FS_ARRAY_SIZE_ERROR: Occurs if ymeas and ymod are not the 
        !!   same length.
        !! - FS_MEMORY_ERROR: Occurs if a memory error is encountered.
    type(single_factor_anova_table) :: rst
        !! A single_factor_anova_table instance containing the ANOVA results.

    ! Local Variables
    integer(int32) :: n, flag
    real(real64), allocatable :: ypack(:)
    real(real64) :: sum_all
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    n = size(ymeas)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    rst%within_factor%f_statistic = ieee_value(sum_all, IEEE_QUIET_NAN)
    rst%within_factor%probability = ieee_value(sum_all, IEEE_QUIET_NAN)

    ! Input Checking
    if (size(ymod) /= n) then
        call report_arrays_not_same_size_error(errmgr, "anova_model_fit", &
            "YMEAS", "YMOD", n, size(ymod))
        return
    end if

    ! Memory Allocation
    allocate(ypack(2 * n), stat = flag)
    if (flag /= 0) then
        call report_memory_error(errmgr, "anova_model_fit", flag)
        return
    end if

    ! Determine the number of DOF
    rst%main_factor%dof = nmodelparams - 1
    rst%within_factor%dof = n - rst%main_factor%dof - 1
    rst%total_dof = n - 1

    ! Process
    ypack(1:n) = ymeas
    ypack(n+1:2*n) = ymod
    rst%overall_mean = mean(ypack)
    rst%total_sum_of_squares = sum((ymeas - rst%overall_mean)**2)
    rst%main_factor%sum_of_squares = sum((ymod - rst%overall_mean)**2)
    rst%within_factor%sum_of_squares = sum((ymeas - ymod)**2)

    rst%total_variance = rst%total_sum_of_squares / rst%total_dof
    rst%main_factor%variance = rst%main_factor%sum_of_squares / &
        rst%main_factor%dof
    rst%within_factor%variance = rst%within_factor%sum_of_squares / &
        rst%within_factor%dof

    ! Compute the F-statistic and probability term
        call anova_probability( &
        rst%main_factor%variance, &
        rst%within_factor%variance, &
        rst%main_factor%dof, &
        rst%within_factor%dof, &
        rst%main_factor%f_statistic, &
        rst%main_factor%probability &
    )    

    ! Formatting
100 format(A, I0, A, I0, A)
101 format(A, I0, A)
end function

! ******************************************************************************
! PRIVATE ROUTINES
! ------------------------------------------------------------------------------
subroutine anova_probability(v1, v2, dof1, dof2, f, p)
    ! Arguments
    real(real64), intent(in) :: v1, v2, dof1, dof2
    real(real64), intent(out) :: f, p

    ! Local Variables
    real(real64) :: d1, d2, a, b, x
    
    ! Process
    f = v1 / v2
    d1 = dof1
    d2 = dof2

    a = 0.5d0 * d2
    b = 0.5d0 * d1
    x = d2 / (d2 + d1 * f)

    p = regularized_beta(a, b, x)
    if (p > 1.0d0) then
        p = 2.0d0 - p
    end if
end subroutine

! ------------------------------------------------------------------------------
end module
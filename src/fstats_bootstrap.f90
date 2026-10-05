! SPDX-FileCopyrightText: 2022-2026 Jason Christopherson
! SPDX-License-Identifier: MIT
!
! MIT License
!
! Copyright (c) 2022-2026 Jason Christopherson
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module fstats_bootstrap
    use iso_fortran_env
    use fstats_errors
    use omp_lib
    use fstats_distributions
    use fstats_descriptive_statistics
    use fstats_special_functions
    use fstats_regression
    use fstats_sampling
    implicit none
    private
    public :: bootstrap_resampling_routine
    public :: bootstrap_statistic_routine
    public :: random_resample
    public :: random_resample_with_replacement
    public :: bootstrap_statistics
    public :: bootstrap

! REFERENCES:
! - https://medium.com/@m21413108/bootstrapping-maximum-entropy-non-parametric-boot-python-3b1e23ea589d
! - https://cran.r-project.org/web/packages/meboot/vignettes/meboot.pdf
! - https://gist.github.com/christianjauregui/314456688a3c2fead43a48be3a47dad6

    type bootstrap_statistics
        !! A collection of statistics resulting from the bootstrap process.
        real(real64) :: statistic_value
            !! The value of the statistic of interest.
        real(real64) :: upper_confidence_interval
            !! The upper confidence limit on the statistic.
        real(real64) :: lower_confidence_interval
            !! The lower confidence limit on the statistic.
        real(real64) :: bias
            !! The bias in the statistic.
        real(real64) :: standard_error
            !! The standard error.
        real(real64), allocatable, dimension(:) :: population
            !! The bootstrap statistic replicates, excluding the observed
            !! statistic.
    end type

    interface
        subroutine bootstrap_resampling_routine(x, xn)
            !! Defines the signature of a subroutine used to compute a 
            !! resampling of data for bootstrapping purposes.
            use iso_fortran_env, only : real64
            real(real64), intent(in), dimension(:) :: x
                !! The N-element array to resample.
            real(real64), intent(out), dimension(size(x)) :: xn
                !! An N-element array where the resampled data set will be 
                !! written.
        end subroutine

        function bootstrap_statistic_routine(x) result(rst)
            !! Defines the signature of a function for computing the desired
            !! bootstrap statistic.
            use iso_fortran_env, only : real64
            real(real64), intent(in), dimension(:) :: x
                !! The array of data to analyze.
            real(real64) :: rst
                !! The resulting statistic.
        end function
    end interface

contains
! ******************************************************************************
! RESAMPLING
! ------------------------------------------------------------------------------
subroutine random_resample(x, xn)
    !! Parametric Gaussian resampling using the sample mean and standard
    !! deviation. This assumes the observations are normally distributed;
    !! it is not the nonparametric bootstrap.
    real(real64), intent(in), dimension(:) :: x
        !! The N-element array to resample.
    real(real64), intent(out), dimension(size(x)) :: xn
        !! An N-element array where the resampled data set will be written.

    ! Local Variables
    integer(int32) :: n
    real(real64) :: avg, sigma

    ! Process
    n = size(x)
    avg = mean(x)
    sigma = standard_deviation(x)
    xn = box_muller_sample(avg, sigma, n)
end subroutine

! ------------------------------------------------------------------------------
subroutine random_resample_with_replacement(x, xn)
    !! Random resampling with replacement from the supplied sample.
    real(real64), intent(in), dimension(:) :: x
        !! The N-element array to resample.
    real(real64), intent(out), dimension(size(x)) :: xn
        !! An N-element array where the resampled data set will be written.

    ! Local Variables
    integer(int32) :: i, n, idx
    real(real64) :: u

    ! Process
    n = size(x)
    do i = 1, n
        call random_number(u)
        idx = int(floor(u * n)) + 1
        if (idx > n) idx = n
        xn(i) = x(idx)
    end do
end subroutine

! ******************************************************************************
! BOOTSTRAPPING
! ------------------------------------------------------------------------------
function bootstrap(stat, x, method, nsamples, alpha) result(rst)
    !! Performs a bootstrap calculation on the supplied data set for the given
    !! statistic.  The default implementation utlizes a random resampling with
    !! replacement.  Other resampling methods may be defined by specifying an 
    !! appropriate routine by means of the method input.  The module provides
    !! another option, [[random_resample]], which is a parametric Gaussian
    !! bootstrap based on the sample mean and standard deviation. It assumes
    !! normally distributed observations. The confidence limits are percentile
    !! intervals calculated from the bootstrap replicates.
    procedure(bootstrap_statistic_routine), pointer, intent(in) :: stat
        !! The routine used to compute the desired statistic.
    real(real64), intent(in), dimension(:) :: x
        !! The N-element data set.
    procedure(bootstrap_resampling_routine), pointer, intent(in), optional :: method
        !! An optional pointer to the method to use for resampling of the data.
        !! If no method is supplied, a random resampling is utilized.
    integer(int32), intent(in), optional :: nsamples
        !! An optional input, that if supplied, specifies the number of
        !! bootstrap replicates to generate. The default is 10 000.
    real(real64), intent(in), optional :: alpha
        !! An optional input, that if supplied, defines the significance level
        !! to use for the analysis.  The default is 0.05.
    type(bootstrap_statistics) :: rst
        !! The resulting bootstrap_statistics type containing the confidence
        !! intervals, bias, standard error, etc. for the analyzed statistic.

    ! Parameters
    real(real64), parameter :: half = 0.5d0
    real(real64), parameter :: p05 = 5.0d-2

    ! Local Variables
    integer(int32) :: i, n, ns
    real(real64) :: a
    real(real64), allocatable, dimension(:) :: xn
    procedure(bootstrap_resampling_routine), pointer :: resample

    ! Initialization
    n = size(x)
    if (present(method)) then
        resample => method
    else
        resample => random_resample_with_replacement
    end if
    if (present(nsamples)) then
        ns = nsamples
    else
        ns = 10000
    end if
    if (present(alpha)) then
        a = alpha
    else
        a = p05
    end if

    if (n < 1 .or. ns < 1) error stop FS_INVALID_INPUT_ERROR
    if (a <= 0.0d0 .or. a >= 1.0d0) error stop FS_INVALID_INPUT_ERROR

    allocate(rst%population(ns))

    ! Analyze the basic data set
    rst%statistic_value = stat(x)

    ! Resampling Process
#ifdef USEOPENMP
    ! Use OpenMP to run operations in parallel
!$OMP PARALLEL DO PRIVATE(xn) SHARED(rst)
    do i = 1, ns
        ! Per-thread memory allocation
        if (.not.allocated(xn)) allocate(xn(n))

        ! Resample the data
        call resample(x, xn)

        ! Compute the statistic
        rst%population(i) = stat(xn)
    end do
!$OMP END PARALLEL DO
#else
    ! OpenMP is not available - run in a serial manner
    allocate(xn(n))
    do i = 1, ns
        ! Resample the data
        call resample(x, xn)

        ! Compute the statistic for the resampled data
        rst%population(i) = stat(xn)
    end do
#endif

    ! Compute the relevant quantities on the resampled statistic
    rst%upper_confidence_interval = quantile(rst%population, 1.0d0 - half * a)
    rst%lower_confidence_interval = quantile(rst%population, half * a)
    rst%bias = mean(rst%population) - rst%statistic_value
    rst%standard_error = standard_deviation(rst%population)
end function

! ------------------------------------------------------------------------------
end module
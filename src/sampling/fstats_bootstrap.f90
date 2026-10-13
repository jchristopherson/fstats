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
    use linalg, only : sort
    implicit none
    private
    public :: bootstrap_resampling_routine
    public :: bootstrap_statistic_routine
    public :: random_resample
    public :: random_resample_with_replacement
    public :: random_resample_paired
    public :: random_resample_clusters
    public :: random_resample_circular_blocks
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

! ------------------------------------------------------------------------------
subroutine random_resample_paired(x, xn)
    !! Resamples paired observations with replacement.
    !! The input stores each pair as adjacent values: [x1, y1, x2, y2, ...].
    real(real64), intent(in), dimension(:) :: x
        !! An even-length array of interleaved paired observations.
    real(real64), intent(out), dimension(size(x)) :: xn
        !! The resampled pairs in the same interleaved layout.

    integer(int32) :: i, n, npairs, idx
    real(real64) :: u

    n = size(x)
    if (n < 2 .or. mod(n, 2) /= 0) error stop FS_INVALID_INPUT_ERROR
    npairs = n / 2
    do i = 1, npairs
        call random_number(u)
        idx = floor(u * npairs, int32) + 1
        idx = min(idx, npairs)
        xn(2 * i - 1:2 * i) = x(2 * idx - 1:2 * idx)
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine random_resample_clusters(x, xn)
    !! Resamples balanced clusters with replacement.
    !! The input stores contiguous records as [cluster_id, value], with every
    !! cluster having the same number of records. The statistic routine should
    !! use the value fields and ignore the cluster identifiers.
    real(real64), intent(in), dimension(:) :: x
        !! An even-length array of contiguous cluster-id/value records.
    real(real64), intent(out), dimension(size(x)) :: xn
        !! The resampled records in the same layout.

    integer(int32) :: i, j, k, n, nrecords, nclusters, cluster_size
    integer(int32) :: source_row, target_row, last_row, idx
    integer(int32), allocatable :: cluster_starts(:)
    real(real64), allocatable :: cluster_ids(:), sorted_ids(:)
    real(real64) :: u

    n = size(x)
    if (n < 2 .or. mod(n, 2) /= 0) error stop FS_INVALID_INPUT_ERROR
    nrecords = n / 2
    allocate(cluster_starts(nrecords), cluster_ids(nrecords), &
        sorted_ids(nrecords))

    nclusters = 1
    cluster_starts(1) = 1
    cluster_ids(1) = x(1)
    do i = 2, nrecords
        if (x(2 * i - 1) /= x(2 * (i - 1) - 1)) then
            nclusters = nclusters + 1
            cluster_starts(nclusters) = i
            cluster_ids(nclusters) = x(2 * i - 1)
        end if
    end do

    cluster_size = nrecords
    if (nclusters > 1) cluster_size = cluster_starts(2) - cluster_starts(1)
    do i = 1, nclusters
        if (i < nclusters) then
            last_row = cluster_starts(i + 1) - 1
        else
            last_row = nrecords
        end if
        if (last_row - cluster_starts(i) + 1 /= cluster_size) &
            error stop FS_INVALID_INPUT_ERROR
    end do

    sorted_ids(1:nclusters) = cluster_ids(1:nclusters)
    call sort(sorted_ids(1:nclusters), .true.)
    do i = 2, nclusters
        if (sorted_ids(i) == sorted_ids(i - 1)) error stop FS_INVALID_INPUT_ERROR
    end do

    do i = 1, nclusters
        call random_number(u)
        idx = floor(u * nclusters, int32) + 1
        idx = min(idx, nclusters)
        source_row = cluster_starts(idx)
        target_row = (i - 1) * cluster_size + 1
        do j = 0, cluster_size - 1
            k = target_row + j
            xn(2 * k - 1:2 * k) = &
                x(2 * (source_row + j) - 1:2 * (source_row + j))
        end do
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine random_resample_circular_blocks(x, xn)
    !! Resamples a time series using a circular moving-block bootstrap.
    !! The block length is floor(sqrt(n)), with a minimum of two when possible.
    !! Define a custom callback to use a domain-specific block length.
    real(real64), intent(in), dimension(:) :: x
        !! The N-element time series.
    real(real64), intent(out), dimension(size(x)) :: xn
        !! The resampled series, assembled from circular blocks.

    integer(int32) :: i, j, n, block_size, nblocks, start, output_index
    integer(int32) :: source_index
    real(real64) :: u

    n = size(x)
    if (n < 1) error stop FS_INVALID_INPUT_ERROR
    block_size = min(n, max(2_int32, int(sqrt(real(n, real64)), int32)))
    nblocks = (n + block_size - 1) / block_size
    output_index = 1
    do i = 1, nblocks
        call random_number(u)
        start = floor(u * n, int32) + 1
        start = min(start, n)
        do j = 0, block_size - 1
            if (output_index > n) exit
            source_index = mod(start - 1 + j, n) + 1
            xn(output_index) = x(source_index)
            output_index = output_index + 1
        end do
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

    ! Compute the relevant quantities on the resampled statistic
    rst%upper_confidence_interval = quantile(rst%population, 1.0d0 - half * a)
    rst%lower_confidence_interval = quantile(rst%population, half * a)
    rst%bias = mean(rst%population) - rst%statistic_value
    rst%standard_error = standard_deviation(rst%population)
end function

! ------------------------------------------------------------------------------
end module
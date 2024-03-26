module fstats_smoothing
    use iso_fortran_env
    use ferror
    use fstats_errors
    use linalg, only : sort
    implicit none
    private
    public :: lowess
contains
! ------------------------------------------------------------------------------
subroutine lowess(x, y, ys, fsmooth, nstps, del, rweights, resid, err)
    !! Computes the smoothing of a data set using a robust locally weighted
    !! scatterplot smoothing (LOWESS) algorithm.  Fitted values are computed at
    !! each of the supplied x values.
    !!
    !! Remarks
    !!
    !! The code is a reimplementation of the LOWESS library.  For a detailed
    !! understanding, see [this]
    !! (http://www.aliquote.org/cours/2012_biomed/biblio/Cleveland1979.pdf) 
    !! paper by William Cleveland.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the independent variable data.  This
        !! array must be monotonically increasing.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the dependent variable data.
    real(real64), intent(out), dimension(:) :: ys
        !! An N-element array where the smoothed results will be written.
    real(real64), intent(in), optional :: fsmooth
        !! An optional input that specifies the amount of smoothing.  
        !! Specifically, this value is the fraction of points used to compute
        !! each value.  As this value increases, the output becomes smoother.
        !! Choosing a value in the range of 0.2 to 0.8 typically results in a
        !! good fit.  The default value is 0.2.
    integer(int32), intent(in), optional :: nstps
        !! An optional input that specifies the numb of iterations.  If set to
        !! zero, a non-robust fit is returned.  The default value is set to 2.
    real(real64), intent(in), optional :: del
        !!
    real(real64), intent(out), optional, dimension(:), target :: rweights
        !! An optional N-element array, that if supplied, will be used to
        !! return the weights given to each data point.
    real(real64), intent(out), optional, dimension(:), target :: resid
        !! An optional N-element array, that if supplied, will be used to 
        !! return the residual.
    class(errors), intent(inout), optional, target :: err
        !! A mechanism for communicating errors and warnings to the 
        !! caller.  Possible warning and error codes are as follows.
        !! - FS_NO_ERROR: No errors encountered.
        !! - FS_ARRAY_SIZE_ERROR: Occurs if any of the arrays are not 
        !!      approriately sized.
        !! - FS_MEMORY_ERROR: Occurs if there is a memory allocation error.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: p2 = 2.0d-1
    real(real64), parameter :: one = 1.0d0
    real(real64), parameter :: three = 3.0d0
    real(real64), parameter :: p001 = 1.0d-3
    real(real64), parameter :: p999 = 0.999d0

    ! Local Variables
    logical :: ok
    integer(int32) :: iter, i, j, nleft, nright, ns, last, m1, m2, n, nsteps, flag
    real(real64) :: f, delta, d1, d2, denom, alpha, cut, eps, cmad, c1, c9, r
    real(real64), allocatable, target, dimension(:) :: rwdef, rsdef
    real(real64), pointer, dimension(:) :: rw, res
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)
    if (present(fsmooth)) then
        f = fsmooth
    else
        f = p2
    end if
    
    if (present(nstps)) then
        nsteps = nstps
    else
        nsteps = 2
    end if

    if (present(del)) then
        delta = del
    else
        delta = 0.0d0
    end if

    if (present(rweights)) then
        if (size(rweights) /= n) then
            call report_array_size_error(errmgr, "lowess", "rweights", n, &
                size(rweights))
            return
        end if
        rw => rweights
    else
        allocate(rwdef(n), stat = flag)
        if (flag /= 0) then
            call report_memory_error(errmgr, "lowess", flag)
            return
        end if
        rw => rwdef
    end if

    if (present(resid)) then
        if (size(resid) /= n) then
            call report_array_size_error(errmgr, "lowess", "resid", n, &
                size(resid))
            return
        end if
        res => resid
    else
        allocate(rsdef(n), stat = flag)
        if (flag /= 0) then
            call report_memory_error(errmgr, "lowess", flag)
            return
        end if
        res => rsdef
    end if
    ns = max(min(int(f * real(n), int32), n), 2)
    eps = epsilon(eps)

    ! Input Checking
    if (size(y) /= n) then
        call report_array_size_error(errmgr, "lowess", "y", n, size(y))
        return
    end if
    if (size(ys) /= n) then
        call report_array_size_error(errmgr, "lowess", "ys", n, size(ys))
        return
    end if

    ! Quick Return
    if (n < 2) then
        ys = y
        return
    end if
    
    ! Process
    do iter = 1, nsteps + 1
        nleft = 1
        nright = ns
        last = 0
        i = 1
        do
            do while (nright < n)
                d1 = x(i) - x(nleft)
                d2 = x(nright+1) - x(i)
                if (d1 <= d2) exit
                nleft = nleft + 1
                nright = nright + 1
            end do

            call lowest(x, y, x(i), ys(i), nleft, nright, res, iter > 1, &
                rw, ok)
            if (.not.ok) ys(i) = y(i)
            if (last < i - 1) then
                denom = x(i) - x(last)
                do j = last + 1, i - 1
                    alpha = (x(j) - x(last)) / denom
                    ys(j) = alpha * ys(i) + (one - alpha) * ys(last)
                end do
            end if
            last = i
            cut = x(last) + delta
            do i = last + 1, n
                if (x(i) > cut) exit
                if (abs(x(i) - x(last)) < eps) then
                    ys(i) = ys(last)
                    last = i
                end if
            end do
            i = max(last + 1, i - 1)

            if (last >= n) exit
        end do

        res = y - ys
        if (iter > nsteps) exit
        rw = abs(res)
        call sort(rw, .true.)
        m1 = 1 + n / 2
        m2 = n - m1 + 1
        cmad = three * (rw(m1) + rw(m2))
        c9 = p999 * cmad
        c1 = p001 * cmad
        do i = 1, n
            r = abs(res(i))
            if (r <= c1) then
                rw(i) = one
            else if (r > c9) then
                rw(i) = zero
            else
                rw(i) = (one - (r / cmad)**2)**2
            end if
        end do
    end do
end subroutine

! ******************************************************************************
! PRIVATE ROUTINES
! ------------------------------------------------------------------------------
! REF:
! - https://en.wikipedia.org/wiki/Local_regression
! - http://www.aliquote.org/cours/2012_biomed/biblio/Cleveland1979.pdf
subroutine lowest(x, y, xs, ys, nleft, nright, w, userw, rw, ok)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x, y, rw ! N ELEMENT
    real(real64), intent(in) :: xs
    real(real64), intent(out) :: ys
    integer(int32), intent(in) :: nleft, nright
    real(real64), intent(out), dimension(:) :: w ! N ELEMENT
    logical, intent(in) :: userw
    logical, intent(out) :: ok

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0
    real(real64), parameter :: p001 = 1.0d-3
    real(real64), parameter :: p999 = 0.999d0

    ! Local Variables
    integer(int32) :: j, n, nrt
    real(real64) :: range, h, h9, h1, a, b, c, r

    ! Initialization
    n = size(x)
    range = x(n) - x(1)
    h = max(xs - x(nleft), x(nright) - xs)
    h9 = p999 * h
    h1 = p001 * h
    a = zero

    ! Process
    do j = nleft, n
        w(j) = zero
        r = abs(x(j) - xs)
        if (r <= h9) then
            if (r > h1) then
                w(j) = (one - (r / h)**3)**3
            else
                w(j) = one
            end if
            if (userw) w(j) = rw(j) * w(j)
            a = a + w(j)
        else if (x(j) > xs) then
            exit
        end if
    end do

    nrt = j - 1
    if (a <= zero) then
        ok = .false.
    else
        ok = .true.
        w(nleft:nrt) = w(nleft:nrt) / a
        if (h > zero) then
            a = zero
            do j = nleft, nrt
                a = a + w(j) * x(j)
            end do
            b = xs - a
            c = zero
            do j = nleft, nrt
                c = c + w(j) * (x(j) - a)**2
            end do
            if (sqrt(c) > p001 * range) then
                b = b / c
                do j = nleft, nrt
                    w(j) = w(j) * (one + b * (x(j) - a))
                end do
            end if
        end if
        ys = zero
        do j = nleft, nrt
            ys = ys + w(j) * y(j)
        end do
    end if
end subroutine

! ------------------------------------------------------------------------------
end module
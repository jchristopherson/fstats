module fstats_msa_tests
    use iso_fortran_env
    use fortran_test_helper
    use fstats
    implicit none

    ! A deterministic study for which every variance component is known in
    ! closed form.  Each cell contains the replicates m - d, m, m + d, where
    ! the cell mean m is exactly additive in the factor effects.  Consequently
    ! every interaction sum of squares vanishes, the within-cell mean square is
    ! exactly d**2, and each main effect variance reduces to the sample
    ! variance of its effects.
    integer(int32), parameter :: nrep = 3
    integer(int32), parameter :: nparts = 5
    integer(int32), parameter :: nops = 3
    integer(int32), parameter :: nlvls = 2
    real(real64), parameter :: dev = 5.0d-1
    real(real64), parameter :: parts(nparts) = [1.0d1, 1.2d1, 1.5d1, 1.1d1, &
        1.7d1]
    real(real64), parameter :: opers(nops) = [0.0d0, 1.0d0, -5.0d-1]
    real(real64), parameter :: lvls(nlvls) = [0.0d0, 2.0d0]

    ! The sample variances of the effects above
    real(real64), parameter :: vpart_exact = 8.5d0
    real(real64), parameter :: voper_exact = 7.0d0 / 1.2d1
    real(real64), parameter :: vlvl_exact = 2.0d0
    real(real64), parameter :: vrepeat_exact = 2.5d-1

contains
! ------------------------------------------------------------------------------
function build_crossed() result(rst)
    real(real64) :: rst(nrep, nparts, nops)

    integer(int32) :: j, k

    do k = 1, nops
        do j = 1, nparts
            rst(1,j,k) = parts(j) + opers(k) - dev
            rst(2,j,k) = parts(j) + opers(k)
            rst(3,j,k) = parts(j) + opers(k) + dev
        end do
    end do
end function

! ------------------------------------------------------------------------------
function build_expanded() result(rst)
    real(real64) :: rst(nrep, nparts, nops, nlvls)

    integer(int32) :: j, k, m

    do m = 1, nlvls
        do k = 1, nops
            do j = 1, nparts
                rst(1,j,k,m) = parts(j) + opers(k) + lvls(m) - dev
                rst(2,j,k,m) = parts(j) + opers(k) + lvls(m)
                rst(3,j,k,m) = parts(j) + opers(k) + lvls(m) + dev
            end do
        end do
    end do
end function

! ------------------------------------------------------------------------------
function test_gauge_rr_crossed() result(rst)
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-10

    ! Local Variables
    real(real64) :: x(nrep, nparts, nops), vgrr, vtotal
    type(gauge_rr_results) :: g

    ! Initialization
    rst = .true.
    x = build_crossed()
    vgrr = vrepeat_exact + voper_exact
    vtotal = vgrr + vpart_exact

    ! An alpha of unity retains the interaction unconditionally
    g = gauge_rr(x, alpha = 1.0d0)

    ! Test 1 - the study type
    if (g%study_type /= GAUGE_RR_CROSSED_STUDY) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_crossed -1"
    end if

    ! Test 2 - the variance components
    if (.not.assert(g%repeatability%variance, vrepeat_exact, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_crossed -2"
    end if
    if (.not.assert(g%operator_by_part%variance, 0.0d0, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_crossed -3"
    end if
    if (.not.assert(g%operator_variation%variance, voper_exact, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_crossed -4"
    end if
    if (.not.assert(g%part_variation%variance, vpart_exact, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_crossed -5"
    end if

    ! Test 3 - the aggregated components
    if (.not.assert(g%reproducibility%variance, voper_exact, tol) .or. &
        .not.assert(g%gauge_rr%variance, vgrr, tol) .or. &
        .not.assert(g%total_variation%variance, vtotal, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_crossed -6"
    end if

    ! Test 4 - the reported percentages
    if (.not.assert(g%gauge_rr%percent_contribution, &
        1.0d2 * vgrr / vtotal, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_crossed -7"
    end if
    if (.not.assert(g%total_variation%percent_study_variation, 1.0d2, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_crossed -8"
    end if
    if (.not.assert(g%repeatability%percent_contribution + &
        g%reproducibility%percent_contribution + &
        g%part_variation%percent_contribution, 1.0d2, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_crossed -9"
    end if

    ! Test 5 - the number of distinct categories
    if (g%distinct_categories /= &
        int(1.41d0 * sqrt(vpart_exact / vgrr), int32)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_crossed -10"
    end if

    ! Test 6 - the terms belonging only to an expanded study are zeroed
    if (.not.assert(g%third_factor_variation%variance, 0.0d0, tol) .or. &
        .not.assert(g%three_way_interaction%variance, 0.0d0, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_crossed -11"
    end if
end function

! ------------------------------------------------------------------------------
function test_gauge_rr_nested() result(rst)
    logical :: rst

    ! In a nested study each operator measures its own parts, so the part term
    ! is the variation of the parts within each operator.  The reference values
    ! are formed directly from the cell means, which the symmetric replicates
    ! leave untouched.

    ! Parameters
    real(real64), parameter :: tol = 1.0d-10

    ! Local Variables
    integer(int32) :: j, k
    real(real64) :: x(nrep, nparts, nops), vpart, vop, cm(nparts, nops), &
        opmean(nops), gm, mspart
    type(gauge_rr_results) :: g

    ! Initialization.  Scaling the operator effect by the part index gives the
    ! operators genuinely different sets of parts.
    rst = .true.
    do k = 1, nops
        do j = 1, nparts
            cm(j,k) = parts(j) + opers(k) * real(j, real64)
            x(1,j,k) = cm(j,k) - dev
            x(2,j,k) = cm(j,k)
            x(3,j,k) = cm(j,k) + dev
        end do
    end do

    gm = sum(cm) / real(nparts * nops, real64)
    do k = 1, nops
        opmean(k) = sum(cm(:,k)) / real(nparts, real64)
    end do
    mspart = 0.0d0
    do k = 1, nops
        do j = 1, nparts
            mspart = mspart + (cm(j,k) - opmean(k))**2
        end do
    end do
    mspart = real(nrep, real64) * mspart / real(nops * (nparts - 1), real64)
    vpart = (mspart - vrepeat_exact) / real(nrep, real64)
    vop = (real(nparts * nrep, real64) * sum((opmean - gm)**2) / &
        real(nops - 1, real64) - mspart) / real(nparts * nrep, real64)

    g = gauge_rr(x, GAUGE_RR_NESTED_STUDY)

    ! Test 1 - the study type, and no interaction is reported
    if (g%study_type /= GAUGE_RR_NESTED_STUDY) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_nested -1"
    end if
    if (g%interaction_pooled .or. &
        .not.assert(g%operator_by_part%variance, 0.0d0, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_nested -2"
    end if

    ! Test 2 - the repeatability is exactly d**2
    if (.not.assert(g%repeatability%variance, vrepeat_exact, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_nested -3"
    end if

    ! Test 3 - the nested part and operator components
    if (.not.assert(g%part_variation%variance, vpart, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_nested -4"
    end if
    if (.not.assert(g%operator_variation%variance, max(vop, 0.0d0), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_nested -5"
    end if

    ! Test 4 - reproducibility is the operator term alone
    if (.not.assert(g%reproducibility%variance, &
        g%operator_variation%variance, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_nested -6"
    end if

    ! Test 5 - the ANOVA degrees of freedom follow the nested design
    if (.not.assert(g%anova_table%operator%dof, &
        real(nops - 1, real64), tol) .or. &
        .not.assert(g%anova_table%part%dof, &
        real(nops * (nparts - 1), real64), tol) .or. &
        .not.assert(g%anova_table%residual%dof, &
        real(nops * nparts * (nrep - 1), real64), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_nested -7"
    end if

    ! Test 6 - the sums of squares decompose the total
    if (.not.assert(g%anova_table%operator%sum_of_squares + &
        g%anova_table%part%sum_of_squares + &
        g%anova_table%residual%sum_of_squares, &
        g%anova_table%total_sum_of_squares, 1.0d-8)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_nested -8"
    end if
end function

! ------------------------------------------------------------------------------
function test_gauge_rr_expanded() result(rst)
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-10

    ! Local Variables
    real(real64) :: x(nrep, nparts, nops, nlvls), vgrr, vtotal, ssum
    type(gauge_rr_results) :: g

    ! Initialization
    rst = .true.
    x = build_expanded()
    vgrr = vrepeat_exact + voper_exact + vlvl_exact
    vtotal = vgrr + vpart_exact

    ! The data are purely additive, so every interaction vanishes and each
    ! main effect variance reduces to the sample variance of its effects
    g = gauge_rr(x, alpha = 1.0d0)

    ! Test 1 - the study type
    if (g%study_type /= GAUGE_RR_EXPANDED_STUDY) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_expanded -1"
    end if

    ! Test 2 - the main effect components
    if (.not.assert(g%repeatability%variance, vrepeat_exact, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_expanded -2"
    end if
    if (.not.assert(g%operator_variation%variance, voper_exact, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_expanded -3"
    end if
    if (.not.assert(g%third_factor_variation%variance, vlvl_exact, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_expanded -4"
    end if
    if (.not.assert(g%part_variation%variance, vpart_exact, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_expanded -5"
    end if

    ! Test 3 - the interactions vanish
    if (.not.assert(g%operator_by_part%variance, 0.0d0, tol) .or. &
        .not.assert(g%third_factor_by_part%variance, 0.0d0, tol) .or. &
        .not.assert(g%third_factor_by_operator%variance, 0.0d0, tol) .or. &
        .not.assert(g%three_way_interaction%variance, 0.0d0, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_expanded -6"
    end if

    ! Test 4 - every term other than the part term is reproducibility
    if (.not.assert(g%reproducibility%variance, &
        voper_exact + vlvl_exact, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_expanded -7"
    end if
    if (.not.assert(g%gauge_rr%variance, vgrr, tol) .or. &
        .not.assert(g%total_variation%variance, vtotal, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_expanded -8"
    end if

    ! Test 5 - the sums of squares decompose the total
    ssum = g%anova_table%part%sum_of_squares + &
        g%anova_table%operator%sum_of_squares + &
        g%anova_table%third_factor%sum_of_squares + &
        g%anova_table%operator_by_part%sum_of_squares + &
        g%anova_table%third_factor_by_part%sum_of_squares + &
        g%anova_table%third_factor_by_operator%sum_of_squares + &
        g%anova_table%three_way%sum_of_squares + &
        g%anova_table%residual%sum_of_squares
    if (.not.assert(ssum, g%anova_table%total_sum_of_squares, 1.0d-8)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_expanded -9"
    end if

    ! Test 6 - the degrees of freedom decompose the total
    if (.not.assert(g%anova_table%part%dof + g%anova_table%operator%dof + &
        g%anova_table%third_factor%dof + &
        g%anova_table%operator_by_part%dof + &
        g%anova_table%third_factor_by_part%dof + &
        g%anova_table%third_factor_by_operator%dof + &
        g%anova_table%three_way%dof + g%anova_table%residual%dof, &
        g%anova_table%total_dof, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_expanded -10"
    end if
end function

! ------------------------------------------------------------------------------
function test_gauge_rr_pooling() result(rst)
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-10

    ! Local Variables
    real(real64) :: x(nrep, nparts, nops), pooled
    type(gauge_rr_results) :: g
    type(gauge_rr_anova_table) :: t

    ! Initialization
    rst = .true.
    x = build_crossed()

    ! The study contains no interaction, so the default behavior must pool the
    ! interaction into the repeatability
    g = gauge_rr(x)
    t = g%anova_table
    pooled = (t%operator_by_part%sum_of_squares + &
        t%residual%sum_of_squares) / (t%operator_by_part%dof + t%residual%dof)

    if (.not.g%interaction_pooled) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_pooling -1"
    end if
    if (.not.assert(g%repeatability%variance, pooled, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_pooling -2"
    end if
    if (.not.assert(g%operator_by_part%variance, 0.0d0, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_pooling -3"
    end if

    ! Retaining the interaction must return the unpooled residual
    g = gauge_rr(x, alpha = 1.0d0)
    if (g%interaction_pooled) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_pooling -4"
    end if
    if (.not.assert(g%repeatability%variance, t%residual%variance, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_pooling -5"
    end if
end function

! ------------------------------------------------------------------------------
function test_gauge_rr_tolerance() result(rst)
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-10
    real(real64), parameter :: trange = 2.0d1
    real(real64), parameter :: k = 5.15d0

    ! Local Variables
    real(real64) :: x(nrep, nparts, nops)
    type(gauge_rr_results) :: g

    ! Initialization
    rst = .true.
    x = build_crossed()

    g = gauge_rr(x, tolerance = trange, multiplier = k, alpha = 1.0d0)

    ! Test 1 - the study multiplier is honored
    if (.not.assert(g%study_multiplier, k, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_tolerance -1"
    end if

    ! Test 2 - the percent tolerance
    if (.not.assert(g%gauge_rr%percent_tolerance, &
        1.0d2 * k * g%gauge_rr%standard_deviation / trange, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_tolerance -2"
    end if

    ! Test 3 - with no tolerance supplied the term is zero
    g = gauge_rr(x, alpha = 1.0d0)
    if (.not.assert(g%gauge_rr%percent_tolerance, 0.0d0, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_gauge_rr_tolerance -3"
    end if
end function

! ------------------------------------------------------------------------------
end module

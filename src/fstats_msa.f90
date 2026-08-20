module fstats_msa
    !! Provides routines supporting measurement systems analysis (MSA).
    use iso_fortran_env
    use fstats_errors
    use fstats_anova, only : anova_factor
    use fstats_special_functions, only : regularized_beta
    implicit none
    private
    public :: GAUGE_RR_CROSSED_STUDY
    public :: GAUGE_RR_NESTED_STUDY
    public :: GAUGE_RR_EXPANDED_STUDY
    public :: gauge_rr_component
    public :: gauge_rr_anova_table
    public :: gauge_rr_results
    public :: gauge_rr

    integer(int32), parameter :: GAUGE_RR_CROSSED_STUDY = 2000
        !! Indicates a crossed study, wherein every operator measures every
        !! part.
    integer(int32), parameter :: GAUGE_RR_NESTED_STUDY = 2001
        !! Indicates a nested study, wherein each operator measures a distinct
        !! set of parts.
    integer(int32), parameter :: GAUGE_RR_EXPANDED_STUDY = 2002
        !! Indicates an expanded study, wherein a third factor is considered
        !! alongside the part and operator factors.

    real(real64), parameter :: default_study_multiplier = 6.0d0
    real(real64), parameter :: default_pooling_alpha = 0.25d0

! ------------------------------------------------------------------------------
    interface gauge_rr
        !! Performs a gauge repeatability and reproducibility (gauge R&R)
        !! study by means of the analysis of variance (ANOVA) method.
        !!
        !! A gauge R&R study apportions the observed variation between the
        !! parts being measured and the measurement system itself.  The
        !! measurement system contribution is split into repeatability, the
        !! variation seen when one operator measures one part repeatedly, and
        !! reproducibility, the variation seen between operators.  Three study
        !! designs are supported.
        !!
        !! ## Crossed
        !!
        !! Every operator measures every part, and each combination is
        !! measured repeatedly.  Because each part is seen by every operator,
        !! the operator-by-part interaction can be separated from the residual,
        !! and the measurements are described by
        !! $$ y_{ijk} = \mu + P_{i} + O_{j} + (PO)_{ij} + \epsilon_{ijk} $$
        !! where \( P \) is the part effect, \( O \) the operator effect,
        !! \( (PO) \) their interaction, and \( \epsilon \) the residual.  The
        !! interaction describes the tendency of operators to differ from one
        !! another by differing amounts from part to part, and it contributes
        !! to reproducibility.  This is the appropriate design whenever the
        !! measurement is non-destructive, and it is the design assumed by
        !! default.  Supply an M-by-P-by-O array and, optionally,
        !! GAUGE_RR_CROSSED_STUDY.
        !!
        !! ## Nested
        !!
        !! Each operator measures its own distinct set of parts, so no part is
        !! shared between operators.  This is the appropriate design when the
        !! measurement destroys or alters the part, or when the parts cannot
        !! otherwise be circulated.  Since no part is common to two operators
        !! there is no information with which to separate an interaction from
        !! the part effect, and the parts are said to be nested within the
        !! operators
        !! $$ y_{ijk} = \mu + O_{j} + P_{i(j)} + \epsilon_{ijk} $$
        !! where \( P_{i(j)} \) denotes the effect of the i-th part measured by
        !! the j-th operator.  Reproducibility is therefore the operator term
        !! alone.  Supply an M-by-P-by-O array together with
        !! GAUGE_RR_NESTED_STUDY; the P parts measured by one operator are
        !! understood to be different from those measured by any other.
        !!
        !! ## Expanded
        !!
        !! A third factor, such as the day, the setup, the fixture, or the
        !! measurement site, is varied alongside the parts and the operators in
        !! a fully crossed arrangement.  All two-way interactions and the
        !! three-way interaction are estimated
        !! $$ y_{ijkl} = \mu + P_{i} + O_{j} + C_{k} + (PO)_{ij} + (PC)_{ik} +
        !! (OC)_{jk} + (POC)_{ijk} + \epsilon_{ijkl} $$
        !! where \( C \) is the third factor.  This design answers whether the
        !! measurement system is stable with respect to an influence that a
        !! conventional study would leave confounded with the repeatability.
        !! Every term other than the part term describes variation introduced
        !! by the measurement system, so all of them contribute to
        !! reproducibility.  Supply an M-by-P-by-O-by-C array; the study type
        !! follows from the rank of the argument.
        module procedure :: gauge_rr_2_factor
        module procedure :: gauge_rr_3_factor
    end interface

! ------------------------------------------------------------------------------
    type gauge_rr_component
        !! Describes a single variance component of a gauge repeatability and
        !! reproducibility study.
        real(real64) :: variance
            !! The estimated variance of the component.
        real(real64) :: standard_deviation
            !! The estimated standard deviation of the component.
        real(real64) :: percent_contribution
            !! The variance of the component expressed as a percentage of the
            !! total variance.  The contributions of the components sum to 100.
        real(real64) :: percent_study_variation
            !! The study variation of the component expressed as a percentage
            !! of the total study variation.  As the study variation is
            !! proportional to a standard deviation rather than a variance,
            !! these percentages do not sum to 100.
        real(real64) :: percent_tolerance
            !! The study variation of the component expressed as a percentage
            !! of the tolerance range.  This term is only meaningful if a
            !! tolerance was supplied; else, it is zero.
    end type

! ------------------------------------------------------------------------------
    type gauge_rr_anova_table
        !! The analysis of variance table underlying a gauge R&R study.  Terms
        !! that do not appear in the study design are zeroed.
        type(anova_factor) :: part
            !! The part term.  In a nested study this term describes the parts
            !! within the operators.
        type(anova_factor) :: operator
            !! The operator term.
        type(anova_factor) :: third_factor
            !! The third factor term of an expanded study.
        type(anova_factor) :: operator_by_part
            !! The operator-by-part interaction.
        type(anova_factor) :: third_factor_by_part
            !! The third-factor-by-part interaction of an expanded study.
        type(anova_factor) :: third_factor_by_operator
            !! The third-factor-by-operator interaction of an expanded study.
        type(anova_factor) :: three_way
            !! The three-way interaction of an expanded study.
        type(anova_factor) :: residual
            !! The residual, or repeatability, term.
        real(real64) :: total_dof
            !! The total number of degrees of freedom.
        real(real64) :: total_sum_of_squares
            !! The total sum of squares.
        real(real64) :: overall_mean
            !! The mean of all of the measurements.
    end type

! ------------------------------------------------------------------------------
    type gauge_rr_results
        !! Contains the results of a gauge repeatability and reproducibility
        !! study.
        integer(int32) :: study_type
            !! The design that was analyzed; one of GAUGE_RR_CROSSED_STUDY,
            !! GAUGE_RR_NESTED_STUDY, or GAUGE_RR_EXPANDED_STUDY.
        type(gauge_rr_component) :: repeatability
            !! The equipment variation, being the variation observed when a
            !! single operator measures a single part repeatedly.
        type(gauge_rr_component) :: reproducibility
            !! The appraiser variation, being the combination of every term
            !! other than the part term and the repeatability.
        type(gauge_rr_component) :: operator_variation
            !! The operator term alone.
        type(gauge_rr_component) :: operator_by_part
            !! The operator-by-part interaction alone.  This term is zero for a
            !! nested study, and zero for a crossed study in which the
            !! interaction was pooled into the repeatability.
        type(gauge_rr_component) :: third_factor_variation
            !! The third factor term alone.  This term is zero unless the study
            !! was expanded.
        type(gauge_rr_component) :: third_factor_by_part
            !! The third-factor-by-part interaction alone.  This term is zero
            !! unless the study was expanded.
        type(gauge_rr_component) :: third_factor_by_operator
            !! The third-factor-by-operator interaction alone.  This term is
            !! zero unless the study was expanded.
        type(gauge_rr_component) :: three_way_interaction
            !! The three-way interaction alone.  This term is zero unless the
            !! study was expanded, and zero if it was pooled into the
            !! repeatability.
        type(gauge_rr_component) :: gauge_rr
            !! The combined repeatability and reproducibility.
        type(gauge_rr_component) :: part_variation
            !! The part-to-part variation.
        type(gauge_rr_component) :: total_variation
            !! The total variation.
        integer(int32) :: distinct_categories
            !! The number of distinct categories the measurement system is able
            !! to resolve.  A value of 5 or more is generally considered
            !! acceptable.
        logical :: interaction_pooled
            !! True if the highest-order interaction was found to be
            !! insignificant and was pooled into the repeatability term.  A
            !! nested study has no interaction, so this term is always false.
        real(real64) :: interaction_probability
            !! The probability associated with the highest-order interaction,
            !! upon which the pooling decision was based.
        real(real64) :: study_multiplier
            !! The multiplier used to convert a standard deviation into a study
            !! variation.
        type(gauge_rr_anova_table) :: anova_table
            !! The ANOVA table from which the variance components were derived.
    end type

contains
! ------------------------------------------------------------------------------
function gauge_rr_2_factor(x, study, tolerance, alpha, multiplier) result(rst)
    !! Performs a crossed or nested gauge repeatability and reproducibility
    !! study.  See the [[gauge_rr]] interface for a description of the designs.
    !!
    !! The variance of each term is estimated from the corresponding mean
    !! squares.  In a crossed study, if the operator-by-part interaction is
    !! found to be insignificant it is pooled into the repeatability term and
    !! the remaining components are re-estimated from the pooled residual.  Any
    !! variance estimate that would otherwise be negative is taken as zero.
    !!
    !! See Also
    !!
    !! - <a href="https://en.wikipedia.org/wiki/ANOVA_gauge_R%26R" target="_blank">Wikipedia</a>
    real(real64), intent(in), dimension(:,:,:) :: x
        !! An M-by-P-by-O array containing the M replicate measurements of each
        !! of the P parts, taken by each of the O operators.  In a nested study
        !! the P parts measured by each operator are distinct from those
        !! measured by any other.  At least two replicates, two parts, and two
        !! operators are required.
    integer(int32), intent(in), optional :: study
        !! An optional input specifying the design, either
        !! GAUGE_RR_CROSSED_STUDY or GAUGE_RR_NESTED_STUDY.  The default is a
        !! crossed study.
    real(real64), intent(in), optional :: tolerance
        !! An optional input specifying the width of the tolerance range of the
        !! characteristic being measured.  If supplied, the study variation of
        !! each component is additionally expressed as a percentage of this
        !! value.
    real(real64), intent(in), optional :: alpha
        !! An optional input specifying the significance level above which the
        !! operator-by-part interaction is considered insignificant and is
        !! pooled into the repeatability.  The default is 0.25.  A value of 1
        !! retains the interaction unconditionally, and a value of 0 pools it
        !! unconditionally.  The term is ignored by a nested study, which has
        !! no interaction.
    real(real64), intent(in), optional :: multiplier
        !! An optional input specifying the number of standard deviations
        !! spanned by the study variation.  The default is 6, which captures
        !! 99.73% of a normally distributed population.
    type(gauge_rr_results) :: rst
        !! The results of the study.

    ! Local Variables
    integer(int32) :: nrep, nparts, nops, design
    real(real64) :: a, k, tol, msres, vrepeat, vop, vpo, vpart

    ! Initialization
    nrep = size(x, 1)
    nparts = size(x, 2)
    nops = size(x, 3)
    design = GAUGE_RR_CROSSED_STUDY
    if (present(study)) design = study
    call resolve_options(alpha, multiplier, tolerance, a, k, tol)

    ! Input Checking
    if (nrep < 2 .or. nparts < 2 .or. nops < 2) then
        error stop FS_UNDERDEFINED_PROBLEM_ERROR
    end if
    if (design /= GAUGE_RR_CROSSED_STUDY .and. &
        design /= GAUGE_RR_NESTED_STUDY) &
    then
        error stop FS_INVALID_INPUT_ERROR
    end if

    ! Build the ANOVA table
    rst%study_type = design
    if (design == GAUGE_RR_NESTED_STUDY) then
        call nested_anova(x, rst%anova_table)
    else
        call crossed_anova(x, rst%anova_table)
    end if

    ! Extract the variance components
    rst%interaction_probability = rst%anova_table%operator_by_part%probability
    rst%interaction_pooled = design == GAUGE_RR_CROSSED_STUDY .and. &
        rst%interaction_probability > a
    if (design == GAUGE_RR_NESTED_STUDY) then
        ! The parts carry the operator effect within them, so the operator
        ! term is referenced to the parts nested inside it
        msres = rst%anova_table%residual%variance
        vrepeat = msres
        vpo = 0.0d0
        vpart = max((rst%anova_table%part%variance - msres) / &
            real(nrep, real64), 0.0d0)
        vop = max((rst%anova_table%operator%variance - &
            rst%anova_table%part%variance) / real(nparts * nrep, real64), &
            0.0d0)
    else if (rst%interaction_pooled) then
        msres = pool(rst%anova_table%operator_by_part, rst%anova_table%residual)
        vrepeat = msres
        vpo = 0.0d0
        vop = max((rst%anova_table%operator%variance - msres) / &
            real(nparts * nrep, real64), 0.0d0)
        vpart = max((rst%anova_table%part%variance - msres) / &
            real(nops * nrep, real64), 0.0d0)
    else
        msres = rst%anova_table%residual%variance
        vrepeat = msres
        vpo = max((rst%anova_table%operator_by_part%variance - msres) / &
            real(nrep, real64), 0.0d0)
        vop = max((rst%anova_table%operator%variance - &
            rst%anova_table%operator_by_part%variance) / &
            real(nparts * nrep, real64), 0.0d0)
        vpart = max((rst%anova_table%part%variance - &
            rst%anova_table%operator_by_part%variance) / &
            real(nops * nrep, real64), 0.0d0)
    end if

    call assemble(rst, vrepeat, vop, vpo, 0.0d0, 0.0d0, 0.0d0, 0.0d0, vpart, &
        k, tol)
end function

! ------------------------------------------------------------------------------
function gauge_rr_3_factor(x, tolerance, alpha, multiplier) result(rst)
    !! Performs an expanded gauge repeatability and reproducibility study, in
    !! which a third factor is crossed with the parts and the operators.  See
    !! the [[gauge_rr]] interface for a description of the design.
    !!
    !! If the three-way interaction is found to be insignificant it is pooled
    !! into the repeatability term, and the remaining components are
    !! re-estimated from the pooled residual.  Any variance estimate that would
    !! otherwise be negative is taken as zero.
    real(real64), intent(in), dimension(:,:,:,:) :: x
        !! An M-by-P-by-O-by-C array containing the M replicate measurements of
        !! each of the P parts, taken by each of the O operators, at each of
        !! the C levels of the third factor.  At least two of each are
        !! required.
    real(real64), intent(in), optional :: tolerance
        !! An optional input specifying the width of the tolerance range of the
        !! characteristic being measured.  If supplied, the study variation of
        !! each component is additionally expressed as a percentage of this
        !! value.
    real(real64), intent(in), optional :: alpha
        !! An optional input specifying the significance level above which the
        !! three-way interaction is considered insignificant and is pooled into
        !! the repeatability.  The default is 0.25.  A value of 1 retains the
        !! interaction unconditionally, and a value of 0 pools it
        !! unconditionally.
    real(real64), intent(in), optional :: multiplier
        !! An optional input specifying the number of standard deviations
        !! spanned by the study variation.  The default is 6, which captures
        !! 99.73% of a normally distributed population.
    type(gauge_rr_results) :: rst
        !! The results of the study.

    ! Local Variables
    integer(int32) :: nrep, nparts, nops, nlvls
    real(real64) :: a, k, tol, msres, ms3, vrepeat, vop, vpo, vc, vpc, voc, &
        v3, vpart

    ! Initialization
    nrep = size(x, 1)
    nparts = size(x, 2)
    nops = size(x, 3)
    nlvls = size(x, 4)
    call resolve_options(alpha, multiplier, tolerance, a, k, tol)

    ! Input Checking
    if (nrep < 2 .or. nparts < 2 .or. nops < 2 .or. nlvls < 2) then
        error stop FS_UNDERDEFINED_PROBLEM_ERROR
    end if

    ! Build the ANOVA table
    rst%study_type = GAUGE_RR_EXPANDED_STUDY
    call expanded_anova(x, rst%anova_table)

    ! Extract the variance components.  The three-way interaction is the
    ! reference term for the two-way interactions unless it was pooled away.
    rst%interaction_probability = rst%anova_table%three_way%probability
    rst%interaction_pooled = rst%interaction_probability > a
    if (rst%interaction_pooled) then
        msres = pool(rst%anova_table%three_way, rst%anova_table%residual)
        ms3 = msres
        v3 = 0.0d0
    else
        msres = rst%anova_table%residual%variance
        ms3 = rst%anova_table%three_way%variance
        v3 = max((ms3 - msres) / real(nrep, real64), 0.0d0)
    end if
    vrepeat = msres
    vpo = max((rst%anova_table%operator_by_part%variance - ms3) / &
        real(nlvls * nrep, real64), 0.0d0)
    vpc = max((rst%anova_table%third_factor_by_part%variance - ms3) / &
        real(nops * nrep, real64), 0.0d0)
    voc = max((rst%anova_table%third_factor_by_operator%variance - ms3) / &
        real(nparts * nrep, real64), 0.0d0)
    vpart = max((rst%anova_table%part%variance - &
        rst%anova_table%operator_by_part%variance - &
        rst%anova_table%third_factor_by_part%variance + ms3) / &
        real(nops * nlvls * nrep, real64), 0.0d0)
    vop = max((rst%anova_table%operator%variance - &
        rst%anova_table%operator_by_part%variance - &
        rst%anova_table%third_factor_by_operator%variance + ms3) / &
        real(nparts * nlvls * nrep, real64), 0.0d0)
    vc = max((rst%anova_table%third_factor%variance - &
        rst%anova_table%third_factor_by_part%variance - &
        rst%anova_table%third_factor_by_operator%variance + ms3) / &
        real(nparts * nops * nrep, real64), 0.0d0)

    call assemble(rst, vrepeat, vop, vpo, vc, vpc, voc, v3, vpart, k, tol)
end function

! ******************************************************************************
! PRIVATE ROUTINES
! ------------------------------------------------------------------------------
subroutine resolve_options(alpha, multiplier, tolerance, a, k, tol)
    !! Applies the default values to the optional study parameters.
    real(real64), intent(in), optional :: alpha, multiplier, tolerance
    real(real64), intent(out) :: a, k, tol

    a = default_pooling_alpha
    k = default_study_multiplier
    tol = 0.0d0
    if (present(alpha)) a = alpha
    if (present(multiplier)) k = multiplier
    if (present(tolerance)) tol = tolerance
    if (k <= 0.0d0) error stop FS_INVALID_INPUT_ERROR
    if (tol < 0.0d0) error stop FS_INVALID_INPUT_ERROR
end subroutine

! ------------------------------------------------------------------------------
pure function pool(a, b) result(rst)
    !! Combines two ANOVA terms into a single mean square.
    type(anova_factor), intent(in) :: a, b
    real(real64) :: rst

    rst = (a%sum_of_squares + b%sum_of_squares) / (a%dof + b%dof)
end function

! ------------------------------------------------------------------------------
pure subroutine set_factor(item, ss, dof, msden, dfden)
    !! Populates an ANOVA term from its sum of squares, and tests it against
    !! the supplied denominator.
    type(anova_factor), intent(out) :: item
    real(real64), intent(in) :: ss
        !! The sum of squares.
    real(real64), intent(in) :: dof
        !! The degrees of freedom.
    real(real64), intent(in) :: msden
        !! The mean square against which to test the term.
    real(real64), intent(in) :: dfden
        !! The degrees of freedom of the denominator.

    item%sum_of_squares = ss
    item%dof = dof
    if (dof > 0.0d0) then
        item%variance = ss / dof
    else
        item%variance = 0.0d0
    end if
    call f_probability(item%variance, msden, dof, dfden, item%f_statistic, &
        item%probability)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine set_null_factor(item)
    !! Zeroes an ANOVA term that does not appear in the study design.
    type(anova_factor), intent(out) :: item

    item%sum_of_squares = 0.0d0
    item%dof = 0.0d0
    item%variance = 0.0d0
    item%f_statistic = 0.0d0
    item%probability = 0.0d0
end subroutine

! ------------------------------------------------------------------------------
pure subroutine f_probability(msnum, msden, dfnum, dfden, f, p)
    !! Computes an F statistic and its upper tail probability.  The tail is
    !! evaluated by way of the reflection I(x; a, b) = 1 - I(1 - x; b, a) so
    !! that small probabilities retain their relative accuracy.
    real(real64), intent(in) :: msnum, msden, dfnum, dfden
    real(real64), intent(out) :: f, p

    if (dfnum <= 0.0d0 .or. dfden <= 0.0d0) then
        f = 0.0d0
        p = 1.0d0
    else if (msden > 0.0d0) then
        f = msnum / msden
        p = regularized_beta(0.5d0 * dfden, 0.5d0 * dfnum, &
            dfden / (dfnum * f + dfden))
    else if (msnum > 0.0d0) then
        f = huge(f)
        p = 0.0d0
    else
        f = 0.0d0
        p = 1.0d0
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine crossed_anova(x, tbl)
    !! Computes the ANOVA table of a crossed study.  The sums of squares are
    !! accumulated about their means rather than from raw sums of squares,
    !! which would be subject to cancellation.
    real(real64), intent(in), dimension(:,:,:) :: x
    type(gauge_rr_anova_table), intent(out) :: tbl

    ! Local Variables
    integer(int32) :: i, j, l, r, p, o
    real(real64) :: gm, ssp, sso, sspo, sse, sst, dfp, dfo, dfpo, dfe, mspo, &
        mse
    real(real64), allocatable :: mp(:), mo(:), cell(:,:)

    ! Initialization
    r = size(x, 1)
    p = size(x, 2)
    o = size(x, 3)
    allocate(mp(p), mo(o), cell(p, o))
    gm = sum(x) / real(size(x), real64)
    do i = 1, p
        mp(i) = sum(x(:,i,:)) / real(r * o, real64)
    end do
    do j = 1, o
        mo(j) = sum(x(:,:,j)) / real(r * p, real64)
        do i = 1, p
            cell(i,j) = sum(x(:,i,j)) / real(r, real64)
        end do
    end do

    ! Sums of squares
    ssp = real(o * r, real64) * sum((mp - gm)**2)
    sso = real(p * r, real64) * sum((mo - gm)**2)
    sspo = 0.0d0
    sse = 0.0d0
    do j = 1, o
        do i = 1, p
            sspo = sspo + (cell(i,j) - mp(i) - mo(j) + gm)**2
            do l = 1, r
                sse = sse + (x(l,i,j) - cell(i,j))**2
            end do
        end do
    end do
    sspo = real(r, real64) * sspo
    sst = sum((x - gm)**2)

    ! Degrees of freedom
    dfp = real(p - 1, real64)
    dfo = real(o - 1, real64)
    dfpo = dfp * dfo
    dfe = real(p * o * (r - 1), real64)

    ! Assemble.  In a random effects model the main effects are referenced to
    ! the interaction, and the interaction to the residual.
    mse = sse / dfe
    mspo = sspo / dfpo
    call set_factor(tbl%part, ssp, dfp, mspo, dfpo)
    call set_factor(tbl%operator, sso, dfo, mspo, dfpo)
    call set_factor(tbl%operator_by_part, sspo, dfpo, mse, dfe)
    call set_factor(tbl%residual, sse, dfe, 0.0d0, 0.0d0)
    call set_null_factor(tbl%third_factor)
    call set_null_factor(tbl%third_factor_by_part)
    call set_null_factor(tbl%third_factor_by_operator)
    call set_null_factor(tbl%three_way)
    tbl%total_sum_of_squares = sst
    tbl%total_dof = real(p * o * r - 1, real64)
    tbl%overall_mean = gm
end subroutine

! ------------------------------------------------------------------------------
subroutine nested_anova(x, tbl)
    !! Computes the ANOVA table of a nested study, in which the parts are
    !! nested within the operators.
    real(real64), intent(in), dimension(:,:,:) :: x
    type(gauge_rr_anova_table), intent(out) :: tbl

    ! Local Variables
    integer(int32) :: i, j, l, r, p, o
    real(real64) :: gm, sso, sspo, sse, sst, dfo, dfpo, dfe, mse, mspo
    real(real64), allocatable :: mo(:), cell(:,:)

    ! Initialization
    r = size(x, 1)
    p = size(x, 2)
    o = size(x, 3)
    allocate(mo(o), cell(p, o))
    gm = sum(x) / real(size(x), real64)
    do j = 1, o
        mo(j) = sum(x(:,:,j)) / real(r * p, real64)
        do i = 1, p
            cell(i,j) = sum(x(:,i,j)) / real(r, real64)
        end do
    end do

    ! Sums of squares
    sso = real(p * r, real64) * sum((mo - gm)**2)
    sspo = 0.0d0
    sse = 0.0d0
    do j = 1, o
        do i = 1, p
            sspo = sspo + (cell(i,j) - mo(j))**2
            do l = 1, r
                sse = sse + (x(l,i,j) - cell(i,j))**2
            end do
        end do
    end do
    sspo = real(r, real64) * sspo
    sst = sum((x - gm)**2)

    ! Degrees of freedom
    dfo = real(o - 1, real64)
    dfpo = real(o * (p - 1), real64)
    dfe = real(p * o * (r - 1), real64)

    ! Assemble.  The nested part term is referenced to the residual, and the
    ! operator term to the parts within it.
    mse = sse / dfe
    mspo = sspo / dfpo
    call set_factor(tbl%operator, sso, dfo, mspo, dfpo)
    call set_factor(tbl%part, sspo, dfpo, mse, dfe)
    call set_factor(tbl%residual, sse, dfe, 0.0d0, 0.0d0)
    call set_null_factor(tbl%operator_by_part)
    call set_null_factor(tbl%third_factor)
    call set_null_factor(tbl%third_factor_by_part)
    call set_null_factor(tbl%third_factor_by_operator)
    call set_null_factor(tbl%three_way)
    tbl%total_sum_of_squares = sst
    tbl%total_dof = real(p * o * r - 1, real64)
    tbl%overall_mean = gm
end subroutine

! ------------------------------------------------------------------------------
subroutine expanded_anova(x, tbl)
    !! Computes the ANOVA table of an expanded study, being a fully crossed
    !! three-factor arrangement of parts, operators, and a third factor.
    real(real64), intent(in), dimension(:,:,:,:) :: x
    type(gauge_rr_anova_table), intent(out) :: tbl

    ! Local Variables
    integer(int32) :: i, j, k, l, r, p, o, c
    real(real64) :: gm, ssp, sso, ssc, sspo, sspc, ssoc, ss3, sse, sst
    real(real64) :: dfp, dfo, dfc, dfpo, dfpc, dfoc, df3, dfe, mse, ms3
    real(real64) :: dnp, dno, dnc, dfnp, dfno, dfnc
    real(real64), allocatable :: mp(:), mo(:), mc(:), mpo(:,:), mpc(:,:), &
        moc(:,:), cell(:,:,:)

    ! Initialization
    r = size(x, 1)
    p = size(x, 2)
    o = size(x, 3)
    c = size(x, 4)
    allocate(mp(p), mo(o), mc(c), mpo(p, o), mpc(p, c), moc(o, c), &
        cell(p, o, c))
    gm = sum(x) / real(size(x), real64)
    do i = 1, p
        mp(i) = sum(x(:,i,:,:)) / real(r * o * c, real64)
    end do
    do j = 1, o
        mo(j) = sum(x(:,:,j,:)) / real(r * p * c, real64)
    end do
    do k = 1, c
        mc(k) = sum(x(:,:,:,k)) / real(r * p * o, real64)
    end do
    do j = 1, o
        do i = 1, p
            mpo(i,j) = sum(x(:,i,j,:)) / real(r * c, real64)
        end do
    end do
    do k = 1, c
        do i = 1, p
            mpc(i,k) = sum(x(:,i,:,k)) / real(r * o, real64)
        end do
        do j = 1, o
            moc(j,k) = sum(x(:,:,j,k)) / real(r * p, real64)
            do i = 1, p
                cell(i,j,k) = sum(x(:,i,j,k)) / real(r, real64)
            end do
        end do
    end do

    ! Main effects
    ssp = real(o * c * r, real64) * sum((mp - gm)**2)
    sso = real(p * c * r, real64) * sum((mo - gm)**2)
    ssc = real(p * o * r, real64) * sum((mc - gm)**2)

    ! Two-way interactions
    sspo = 0.0d0
    do j = 1, o
        do i = 1, p
            sspo = sspo + (mpo(i,j) - mp(i) - mo(j) + gm)**2
        end do
    end do
    sspo = real(c * r, real64) * sspo

    sspc = 0.0d0
    ssoc = 0.0d0
    do k = 1, c
        do i = 1, p
            sspc = sspc + (mpc(i,k) - mp(i) - mc(k) + gm)**2
        end do
        do j = 1, o
            ssoc = ssoc + (moc(j,k) - mo(j) - mc(k) + gm)**2
        end do
    end do
    sspc = real(o * r, real64) * sspc
    ssoc = real(p * r, real64) * ssoc

    ! The three-way interaction and the residual
    ss3 = 0.0d0
    sse = 0.0d0
    do k = 1, c
        do j = 1, o
            do i = 1, p
                ss3 = ss3 + (cell(i,j,k) - mpo(i,j) - mpc(i,k) - moc(j,k) + &
                    mp(i) + mo(j) + mc(k) - gm)**2
                do l = 1, r
                    sse = sse + (x(l,i,j,k) - cell(i,j,k))**2
                end do
            end do
        end do
    end do
    ss3 = real(r, real64) * ss3
    sst = sum((x - gm)**2)

    ! Degrees of freedom
    dfp = real(p - 1, real64)
    dfo = real(o - 1, real64)
    dfc = real(c - 1, real64)
    dfpo = dfp * dfo
    dfpc = dfp * dfc
    dfoc = dfo * dfc
    df3 = dfp * dfo * dfc
    dfe = real(p * o * c * (r - 1), real64)
    mse = sse / dfe
    ms3 = ss3 / df3

    ! In a fully random three-factor model no single mean square supplies the
    ! expected value of a main effect, so the denominators are synthesized
    ! from the relevant two-way terms and their degrees of freedom follow from
    ! the Satterthwaite approximation.
    call satterthwaite(sspo / dfpo, dfpo, sspc / dfpc, dfpc, ms3, df3, dnp, &
        dfnp)
    call satterthwaite(sspo / dfpo, dfpo, ssoc / dfoc, dfoc, ms3, df3, dno, &
        dfno)
    call satterthwaite(sspc / dfpc, dfpc, ssoc / dfoc, dfoc, ms3, df3, dnc, &
        dfnc)

    call set_factor(tbl%part, ssp, dfp, dnp, dfnp)
    call set_factor(tbl%operator, sso, dfo, dno, dfno)
    call set_factor(tbl%third_factor, ssc, dfc, dnc, dfnc)
    call set_factor(tbl%operator_by_part, sspo, dfpo, ms3, df3)
    call set_factor(tbl%third_factor_by_part, sspc, dfpc, ms3, df3)
    call set_factor(tbl%third_factor_by_operator, ssoc, dfoc, ms3, df3)
    call set_factor(tbl%three_way, ss3, df3, mse, dfe)
    call set_factor(tbl%residual, sse, dfe, 0.0d0, 0.0d0)
    tbl%total_sum_of_squares = sst
    tbl%total_dof = real(p * o * c * r - 1, real64)
    tbl%overall_mean = gm
end subroutine

! ------------------------------------------------------------------------------
pure subroutine satterthwaite(ms1, df1, ms2, df2, ms3, df3, ms, df)
    !! Forms the synthetic denominator MS1 + MS2 - MS3 along with its effective
    !! degrees of freedom by way of the Satterthwaite approximation.
    real(real64), intent(in) :: ms1, df1, ms2, df2, ms3, df3
    real(real64), intent(out) :: ms, df

    ! Local Variables
    real(real64) :: denom

    ms = ms1 + ms2 - ms3
    denom = 0.0d0
    if (df1 > 0.0d0) denom = denom + ms1**2 / df1
    if (df2 > 0.0d0) denom = denom + ms2**2 / df2
    if (df3 > 0.0d0) denom = denom + ms3**2 / df3
    if (ms > 0.0d0 .and. denom > 0.0d0) then
        df = ms**2 / denom
    else
        ms = 0.0d0
        df = 0.0d0
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine assemble(rst, vrepeat, vop, vpo, vc, vpc, voc, v3, vpart, k, &
    tol)
    !! Aggregates the variance components and populates the reported terms.
    type(gauge_rr_results), intent(inout) :: rst
    real(real64), intent(in) :: vrepeat, vop, vpo, vc, vpc, voc, v3, vpart, &
        k, tol

    ! Local Variables
    real(real64) :: vreprod, vgrr, vtotal

    ! Every term other than the part term describes variation introduced by
    ! the measurement system
    vreprod = vop + vpo + vc + vpc + voc + v3
    vgrr = vrepeat + vreprod
    vtotal = vgrr + vpart

    rst%study_multiplier = k
    call set_component(rst%repeatability, vrepeat, vtotal, k, tol)
    call set_component(rst%reproducibility, vreprod, vtotal, k, tol)
    call set_component(rst%operator_variation, vop, vtotal, k, tol)
    call set_component(rst%operator_by_part, vpo, vtotal, k, tol)
    call set_component(rst%third_factor_variation, vc, vtotal, k, tol)
    call set_component(rst%third_factor_by_part, vpc, vtotal, k, tol)
    call set_component(rst%third_factor_by_operator, voc, vtotal, k, tol)
    call set_component(rst%three_way_interaction, v3, vtotal, k, tol)
    call set_component(rst%gauge_rr, vgrr, vtotal, k, tol)
    call set_component(rst%part_variation, vpart, vtotal, k, tol)
    call set_component(rst%total_variation, vtotal, vtotal, k, tol)

    ! The number of distinct categories is the number of non-overlapping
    ! confidence intervals of the measurement system that span the part
    ! variation, truncated to an integer.  The coefficient is the rounded
    ! value of sqrt(2) published by the AIAG.
    if (vgrr <= 0.0d0) then
        rst%distinct_categories = huge(rst%distinct_categories)
    else
        rst%distinct_categories = int(1.41d0 * sqrt(vpart / vgrr), int32)
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine set_component(item, v, vtotal, k, tol)
    !! Populates a variance component from its variance.
    type(gauge_rr_component), intent(out) :: item
    real(real64), intent(in) :: v
        !! The variance of the component.
    real(real64), intent(in) :: vtotal
        !! The total variance.
    real(real64), intent(in) :: k
        !! The study variation multiplier.
    real(real64), intent(in) :: tol
        !! The tolerance range, or zero if none was supplied.

    item%variance = v
    item%standard_deviation = sqrt(v)
    if (vtotal > 0.0d0) then
        item%percent_contribution = 1.0d2 * v / vtotal
        item%percent_study_variation = 1.0d2 * sqrt(v / vtotal)
    else
        item%percent_contribution = 0.0d0
        item%percent_study_variation = 0.0d0
    end if
    if (tol > 0.0d0) then
        item%percent_tolerance = 1.0d2 * k * item%standard_deviation / tol
    else
        item%percent_tolerance = 0.0d0
    end if
end subroutine

! ------------------------------------------------------------------------------
end module
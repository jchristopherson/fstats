module fstats_statistics_tests
    use iso_fortran_env
    use fstats
    use fstats_test_helper
    implicit none
contains
! ------------------------------------------------------------------------------
    function mean_test_1() result(rst)
        ! Arguments
        logical :: rst

        ! Variables
        integer(int32), parameter :: n = 50
        real(real64) :: x(n), avg, ans

        ! Initialization
        rst = .true.
        x = create_real64_array(n)

        ! Determine the solution
        ans = sum(x) / n

        ! Test 1
        avg = mean(x)
        if (.not.is_equal(avg, ans)) then
            rst = .false.
            print '(A)', "TEST FAILED: Mean Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function variance_test_1() result(rst)
        ! Arguments
        logical :: rst

        ! Variables
        integer(int32), parameter :: n = 50
        real(real64) :: x(n), v, ans, avg

        ! Initialization
        rst = .true.
        x = create_real64_array(n)

        avg = mean(x)

        ans = sum((x - avg)**2) / (n - 1.0d0)

        ! Test 1
        v = variance(x)
        if (.not.is_equal(v, ans)) then
            rst = .false.
            print '(A)', "TEST FAILED: Variance Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function standard_deviation_test_1() result(rst)
        ! Arguments
        logical :: rst

        ! Variables
        integer(int32), parameter :: n = 50
        real(real64) :: x(n), v, ans, avg

        ! Initialization
        rst = .true.
        x = create_real64_array(n)

        avg = mean(x)

        ans = sqrt(sum((x - avg)**2) / (n - 1.0d0))

        ! Test 1
        v = standard_deviation(x)
        if (.not.is_equal(v, ans)) then
            rst = .false.
            print '(A)', "TEST FAILED: Standard Deviation Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function median_test_1() result(rst)
        ! Arguments
        logical :: rst

        ! Variables
        integer(int32), parameter :: n = 7
        real(real64) :: x(n), med, ans

        ! Initialization
        rst = .true.
        x = [1.0d0, 2.0d0, 2.0d0, 4.0d0, 7.0d0, 9.0d0, 3.0d0]
        ans = 3.0d0

        ! Test 1
        med = median(x)
        if (.not.is_equal(ans, med)) then
            rst = .false.
            print '(A)', "TEST FAILED: Median Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function r_squared_test_1() result(rst)
        ! Arguments
        logical :: rst

        ! NOTE:
        ! Data and results computed using Excel

        ! Variables
        integer(int32), parameter :: n = 41
        integer(int32), parameter :: p = 2
        real(real64), parameter :: r2ans = 0.99964955179438d0
        real(real64), parameter :: r2adjans = 0.999631107152d0
        real(real64) :: x(n), xm(n), r2, r2adj

        ! Initialization
        rst = .true.
        x = [ &
            0.74293343964928d0, &
            1.03824302667136d0, &
            1.24885811037080d0, &
            1.41310201321365d0, &
            1.79883628670055d0, &
            2.19087133472520d0, &
            2.31525552428746d0, &
            2.68107856198154d0, &
            2.98235140522939d0, &
            3.41747699679949d0, &
            3.70273035803437d0, &
            3.87952645198981d0, &
            4.19213205695136d0, &
            4.56471706308143d0, &
            4.93604973694783d0, &
            5.08488185540248d0, &
            5.53125594637151d0, &
            5.62434266260204d0, &
            6.09355667145189d0, &
            6.27253567086689d0, &
            6.56795648403612d0, &
            7.03277416903707d0, &
            7.31963808923275d0, &
            7.62307787707421d0, &
            7.80658758580811d0, &
            8.11155981519294d0, &
            8.42961829875520d0, &
            8.72313187683920d0, &
            9.02401354663279d0, &
            9.34923432013754d0, &
            9.59417252155719d0, &
            9.96313533736369d0, &
            10.14422280277640d0, &
            10.52924827792010d0, &
            10.78865600915350d0, &
            11.17866543670860d0, &
            11.39583411689620d0, &
            11.76097700380790d0, &
            11.95842221979940d0, &
            12.39828762076430d0, &
            12.56757580387910d0 &
        ]
        xm = [ &
            0.646498523110604d0, &
            0.945853504743733d0, &
            1.245208486376860d0, &
            1.544563468009990d0, &
            1.843918449643120d0, &
            2.143273431276250d0, &
            2.442628412909380d0, &
            2.741983394542510d0, &
            3.041338376175640d0, &
            3.340693357808760d0, &
            3.640048339441890d0, &
            3.939403321075020d0, &
            4.238758302708150d0, &
            4.538113284341280d0, &
            4.837468265974410d0, &
            5.136823247607540d0, &
            5.436178229240670d0, &
            5.735533210873800d0, &
            6.034888192506930d0, &
            6.334243174140060d0, &
            6.633598155773180d0, &
            6.932953137406310d0, &
            7.232308119039440d0, &
            7.531663100672570d0, &
            7.831018082305700d0, &
            8.130373063938830d0, &
            8.429728045571960d0, &
            8.729083027205090d0, &
            9.028438008838220d0, &
            9.327792990471350d0, &
            9.627147972104480d0, &
            9.926502953737610d0, &
            10.225857935370700d0, &
            10.525212917003900d0, &
            10.824567898637000d0, &
            11.123922880270100d0, &
            11.423277861903300d0, &
            11.722632843536400d0, &
            12.021987825169500d0, &
            12.321342806802600d0, &
            12.620697788435800d0 &    
        ]

        ! Test 1
        r2 = r_squared(x, xm)
        if (.not.is_equal(r2, r2ans)) then
            rst = .false.
            print '(A)', "TEST FAILED: R-Squared Test 1"
        end if
        
        ! Test 2
        r2adj = adjusted_r_squared(p, x, xm)
        if (.not.is_equal(r2adj, r2adjans)) then
            rst = .false.
            print '(A)', "TEST FAILED: R-Squared Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function t_test_test_1() result(rst)
        ! Arguments
        logical :: rst

        ! NOTE:
        ! Data and results computed using Excel

        ! Variables
        integer(int32), parameter :: n1 = 24
        integer(int32), parameter :: n2 = 30
        real(real64), parameter :: p1ans = 0.400311271833311d0
        real(real64), parameter :: p2ans = 0.393261072269881d0
        real(real64), parameter :: p3ans = 0.336451697935195d0
        real(real64) :: x1(n1), x2(n2), t1, p1, d1, t2, p2, d2, t3, p3, d3

        ! Initialization
        rst = .true.
        x1 = [ &
            0.237002480020531d0, &
            0.249842417956438d0, &
            0.251423411811850d0, &
            0.382917286246487d0, &
            0.528741553829583d0, &
            0.457033804084313d0, &
            0.353101298495923d0, &
            0.992831323993046d0, &
            0.665780731344508d0, &
            0.862169232875625d0, &
            0.030392104252987d0, &
            0.782793491156521d0, &
            0.681542313374514d0, &
            0.740264397548036d0, &
            0.349663412525814d0, &
            0.883557553697622d0, &
            0.281040873254668d0, &
            0.519922301130437d0, &
            0.693781206154896d0, &
            0.233070564087244d0, &
            0.599018512926571d0, &
            0.643730371805621d0, &
            0.387975220333301d0, &
            0.030592274330406d0 &
        ]
        x2 = [ &
            0.994458762426464d0, &
            0.166766332238246d0, &
            0.751851412717990d0, &
            0.888805979260137d0, &
            0.372638192453291d0, &
            0.525328139741015d0, &
            0.182360758006977d0, &
            0.405844035648620d0, &
            0.621708871220267d0, &
            0.649565811259767d0, &
            0.802542182211269d0, &
            0.616035912118075d0, &
            0.710391627657817d0, &
            0.943094087018965d0, &
            0.213641637441068d0, &
            0.506435408965975d0, &
            0.108477598968561d0, &
            0.418708813069064d0, &
            0.831979924932566d0, &
            0.603398662772554d0, &
            0.873933396679142d0, &
            0.006005935739189d0, &
            0.688774678554933d0, &
            0.887322822903849d0, &
            0.782795872250114d0, &
            0.175820587186482d0, &
            0.956521100095116d0, &
            0.816615602416826d0, &
            0.226929745400446d0, &
            0.066444929423269d0 &    
        ]

        ! Test 1 - Equal Variances
        call t_test_equal_variance(x1, x2, t1, p1, d1)
        if (.not.is_equal(p1, p1ans)) then
            rst = .false.
            print '(A)', "TEST FAILED: T-Test 1"
        end if

        ! Test 2 - Unequal Variances
        call t_test_unequal_variance(x1, x2, t2, p2, d2)
        if (.not.is_equal(p2, p2ans)) then
            rst = .false.
            print '(A)', "TEST FAILED: T-Test 2"
        end if

        ! Test 3 - Paired
        call t_test_paired(x1, x2(1:n1), t3, p3, d3)
        if (.not.is_equal(p3, p3ans)) then
            rst = .false.
            print '(A)', "TEST FAILED: T-Test 3"
        end if
    end function

! ------------------------------------------------------------------------------
    function f_test_test_1() result(rst)
        ! Arguments
        logical :: rst

        ! Variables
        integer(int32), parameter :: n1 = 25
        integer(int32), parameter :: n2 = 20
        real(real64), parameter :: p1ans = 0.917845200397729d0
        real(real64) :: x1(n1), x2(n2), f1, p1, d1a, d1b

        ! Initialization
        rst = .true.
        x1 = [ &
            0.628237333815969d0, &
            0.880814843302510d0, &
            0.440670600863931d0, &
            0.910375448266385d0, &
            0.561689246761569d0, &
            0.229990505241836d0, &
            0.841676467491122d0, &
            0.677924722574395d0, &
            0.442739209907264d0, &
            0.459747728939493d0, &
            0.168218225021278d0, &
            0.047240720722516d0, &
            0.361664920289616d0, &
            0.680756309597675d0, &
            0.231433224697163d0, &
            0.356021046801774d0, &
            0.168067327623295d0, &
            0.399732643362505d0, &
            0.790399722839469d0, &
            0.132406412559423d0, &
            0.383477756184866d0, &
            0.948900448809131d0, &
            0.125709913245416d0, &
            0.250996116991792d0, &
            0.483364480129666d0 &    
        ]
        x2 = [ &
            0.706493956662183d0, &
            0.450534624824412d0, &
            0.646236078478981d0, &
            0.293181230330176d0, &
            0.750623405641998d0, &
            0.060705272933877d0, &
            0.128825317596501d0, &
            0.449871825416915d0, &
            0.352769551621191d0, &
            0.787163388992404d0, &
            0.271577140729617d0, &
            0.595121914652586d0, &
            0.738831480277178d0, &
            0.368106521227665d0, &
            0.185871577824769d0, &
            0.907905144429555d0, &
            0.313952867238934d0, &
            0.998254961034461d0, &
            0.428081998768432d0, &
            0.113662531842427d0 &    
        ]

        ! Test 1
        call f_test(x1, x2, f1, p1, d1a, d1b)
        if (.not.is_equal(p1, p1ans)) then
            rst = .false.
            print '(A)', "TEST FAILED: F Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function anova_test_1() result(rst)
        ! Arguments
        logical :: rst

        ! This test is from data I measured during a shooting experiment.  The 
        ! reference values were computed using JASP.

        ! Variables
        integer(int32), parameter :: n = 10
        real(real64), parameter :: tol = 1.0d-3
        real(real64), parameter :: bvarans = 352.8d0
        real(real64), parameter :: ivarans = 164.344d0
        real(real64), parameter :: fans = 2.147d0
        real(real64), parameter :: pans = 0.160d0
        real(real64) :: x(n, 2)
        type(single_factor_anova_table) :: rst1

        ! Initialization
        rst = .true.
        x = reshape( &
            [ &
                3.086d3, 3.082d3, 3.069d3, 3.072d3, 3.045d3, 3.070d3, 3.079d3, &
                3.050d3, 3.062d3, 3.062d3, 3.075d3, 3.061d3, 3.063d3, 3.038d3, &
                3.070d3, 3.062d3, 3.070d3, 3.049d3, 3.042d3, 3.063d3 &
            ], &
            [n, 2] &
        )

        ! Test 1
        rst1 = anova(x)
        if (.not.is_equal(rst1%main_factor%variance, bvarans, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: ANOVA 1 - 1"
        end if
        if (.not.is_equal(rst1%within_factor%variance, ivarans, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: ANOVA 1 - 2"
        end if
        if (.not.is_equal(rst1%main_factor%f_statistic, fans, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: ANOVA 1 - 3"
        end if
        if (.not.is_equal(rst1%main_factor%probability, pans, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: ANOVA 1 - 4"
        end if
    end function

! ------------------------------------------------------------------------------
    function anova_test_2() result(rst)
        ! Arguments
        logical :: rst

        ! This test is from data I measured during a shooting experiment.  The 
        ! reference values were computed using JASP.

        ! Variables
        integer(int32), parameter :: ngroups = 2
        integer(int32), parameter :: nshots = 5
        integer(int32), parameter :: nshooters = 2
        real(real64), parameter :: tol = 1.0d-3
        real(real64), parameter :: mf1var = 352.8d0
        real(real64), parameter :: mf1f = 2.003d0
        real(real64), parameter :: mf1p = 0.176d0
        real(real64), parameter :: mf2var = 135.20d0
        real(real64), parameter :: mf2f = 0.768d0
        real(real64), parameter :: mf2p = 0.394d0
        real(real64), parameter :: intvar = 5.0d0
        real(real64), parameter :: intf = 0.028d0
        real(real64), parameter :: intp = 0.868d0
        real(real64), parameter :: evar = 176.125d0
        real(real64) :: x(nshots, ngroups, nshooters)
        type(two_factor_anova_table) :: rst1

        ! Initialization
        rst = .true.
        x(:,:,1) = reshape( &
            [ &
                3.086d3, 3.082d3, 3.069d3, 3.072d3, 3.045d3, 3.075d3, 3.061d3, &
                3.063d3, 3.038d3, 3.070d3 &
            ], &
            [nshots, ngroups] &
        )
        x(:,:,2) = reshape( &
            [ &
                3.070d3, 3.079d3, 3.050d3,  3.062d3, 3.062d3, 3.062d3, 3.070d3, &
                3.049d3, 3.042d3, 3.063d3 &
            ], &
            [nshots, ngroups] &
        )

        ! Test 1
        rst1 = anova(x)
        if (.not.is_equal(rst1%main_factor_1%variance, mf1var, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: ANOVA 2 - 1"
        end if
        if (.not.is_equal(rst1%main_factor_1%f_statistic, mf1f, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: ANOVA 2 - 2"
        end if
        if (.not.is_equal(rst1%main_factor_1%probability, mf1p, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: ANOVA 2 - 3"
        end if
        if (.not.is_equal(rst1%main_factor_2%variance, mf2var, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: ANOVA 2 - 4"
        end if
        if (.not.is_equal(rst1%main_factor_2%f_statistic, mf2f, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: ANOVA 2 - 5"
        end if
        if (.not.is_equal(rst1%main_factor_2%probability, mf2p, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: ANOVA 2 - 6"
        end if
        if (.not.is_equal(rst1%interaction%variance, intvar, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: ANOVA 2 - 7"
        end if
        if (.not.is_equal(rst1%interaction%f_statistic, intf, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: ANOVA 2 - 8"
        end if
        if (.not.is_equal(rst1%interaction%probability, intp, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: ANOVA 2 - 9"
        end if
        if (.not.is_equal(rst1%within_factor%variance, evar, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: ANOVA 2 - 10"
        end if
    end function

! ------------------------------------------------------------------------------
    function anova_test_3() result(rst)
        ! Arguments
        logical :: rst

        ! The reference results were computed using JASP.

        ! Variables
        integer(int32), parameter :: n = 31
        real(real64), parameter :: tol = 1.0d-3
        real(real64), parameter :: mvar = 5.871d0
        real(real64), parameter :: mf = 10226.742d0
        real(real64), parameter :: wssq = 0.017d0
        real(real64), parameter :: tssq = 5.888d0
        real(real64) :: ymeas(n), ymod(n)
        type(single_factor_anova_table) :: tbl

        ! Initialization
        rst = .true.
        ymeas = [ &
            0.57785513826045d0, &
            0.61488309560422d0, &
            0.63389112748856d0, &
            0.71840582970172d0, &
            0.75366850275911d0, &
            0.81496785731015d0, &
            0.86187099649970d0, &
            0.92510053374438d0, &
            0.94703801852006d0, &
            1.02519804334328d0, &
            1.04214235449761d0, &
            1.12152856678444d0, &
            1.17757031499407d0, &
            1.22923756752537d0, &
            1.26111406259387d0, &
            1.29640816255143d0, &
            1.39435365705112d0, &
            1.36714439156037d0, &
            1.42816443143515d0, &
            1.54894493507327d0, &
            1.50510014928299d0, &
            1.56070102375152d0, &
            1.60911301248153d0, &
            1.66368736687550d0, &
            1.70714954545687d0, &
            1.80093594761811d0, &
            1.84781998890644d0, &
            1.88424282167581d0, &
            1.96617423937314d0, &
            1.97700526644311d0, &
            2.03413725715414d0 &
        ]
        ymod = [ &
            0.57000952115854d0, &
            0.61866599223468d0, &
            0.66732246331081d0, &
            0.71597893438695d0, &
            0.76463540546308d0, &
            0.81329187653922d0, &
            0.86194834761535d0, &
            0.91060481869149d0, &
            0.95926128976762d0, &
            1.00791776084375d0, &
            1.05657423191989d0, &
            1.10523070299602d0, &
            1.15388717407216d0, &
            1.20254364514829d0, &
            1.25120011622443d0, &
            1.29985658730056d0, &
            1.34851305837670d0, &
            1.39716952945283d0, &
            1.44582600052897d0, &
            1.49448247160510d0, &
            1.54313894268124d0, &
            1.59179541375737d0, &
            1.64045188483351d0, &
            1.68910835590964d0, &
            1.73776482698578d0, &
            1.78642129806191d0, &
            1.83507776913805d0, &
            1.88373424021418d0, &
            1.93239071129032d0, &
            1.98104718236645d0, &
            2.02970365344259d0 &
        ]

        ! Test 1
        tbl = anova(2, ymeas, ymod)

        if (.not.is_equal(tbl%main_factor%variance, mvar, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: ANOVA 3 - 1"
        end if
        if (.not.is_equal(tbl%main_factor%f_statistic, mf, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: ANOVA 3 - 2"
        end if
        if (.not.is_equal(tbl%within_factor%sum_of_squares, wssq, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: ANOVA 3 - 4"
        end if
        if (.not.is_equal(tbl%total_sum_of_squares, tssq, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: ANOVA 3 - 4"
        end if
    end function

! ------------------------------------------------------------------------------
    function confidence_interval_test_1() result(rst)
        ! Arguments
        logical :: rst

        ! Variables
        real(real64), parameter :: alpha = 2.0d0
        real(real64), parameter :: tol = 1.0d-3
        real(real64), parameter :: t20 = 2.086d0
        real(real64), parameter :: t10 = 2.228d0
        real(real64), parameter :: t30 = 2.042d0
        real(real64), parameter :: sigma = 1.0d0
        integer(int32), parameter :: n = 1
        type(t_distribution) :: dist
        integer(int32) :: flag
        real(real64) :: c

        ! Initialization
        rst = .true.

        ! Test 1, nu = 10
        dist%dof = 10.0d0
        c = confidence_interval(dist, alpha, sigma, n)
        if (.not.is_equal(c, t10, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: Confidence Interval Test 1"
        end if

        ! Test 2, nu = 20
        dist%dof = 20.0d0
        c = confidence_interval(dist, alpha, sigma, n)
        if (.not.is_equal(c, t20, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: Confidence Interval Test 2"
        end if

        ! Test 3, nu = 30
        dist%dof = 30.0d0
        c = confidence_interval(dist, alpha, sigma, n)
        if (.not.is_equal(c, t30, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: Confidence Interval Test 3"
        end if
    end function

! ------------------------------------------------------------------------------
    function beta_test_1() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: a1 = 2.0d0
        real(real64), parameter :: b1 = 2.0d0
        real(real64), parameter :: a3 = 4.0d0
        real(real64), parameter :: b3 = 3.0d0
        real(real64), parameter :: a5 = 1.0d0
        real(real64), parameter :: b5 = 1.0d0

        ! Local Variables
        real(real64) :: dans

        ! Initialization
        rst = .true.

        ! Test 1: a = 2, b = 2: beta = 1/6
        dans = beta(a1, b1)
        if (.not.is_equal(dans, 1.0d0 / 6.0d0)) then
            rst = .false.
            print '(A)', "TEST FAILED: Beta Test 1 - 1"
        end if

        ! Test 2: beta = 1/60
        dans = beta(a3, b3)
        if (.not.is_equal(dans, 1.0d0 / 6.0d1)) then
            rst = .false.
            print '(A)', "TEST FAILED: Beta Test 1 - 2"
        end if

        ! Test 3: beta = 1
        dans = beta(a5, b5)
        if (.not.is_equal(dans, 1.0d0)) then
            rst = .false.
            print '(A)', "TEST FAILED: Beta Test 1 - 3"
        end if
    end function

! ------------------------------------------------------------------------------
    function beta_test_2() result(rst)
        ! Arguments
        logical :: rst

        ! REF: https://keisan.casio.com/exec/system/1180573396
    
        ! Parameters
        real(real64), parameter :: a1 = 1.0d0
        real(real64), parameter :: b1 = 3.0d0
        real(real64), parameter :: x1 = 0.4d0
        real(real64), parameter :: ans1 = 0.2613333333333333333333d0
        real(real64), parameter :: ans1r = 0.784d0

        ! Local Variables
        real(real64) :: dans, dansr

        ! Initialization
        rst = .true.

        ! Test 1
        dans = incomplete_beta(a1, b1, x1)
        if (.not.is_equal(dans, ans1)) then
            rst = .false.
            print '(A)', "TEST FAILED: Beta Test 2 - 1"
        end if

        ! Test 3
        dansr = regularized_beta(a1, b1, x1)
        if (.not.is_equal(dansr, ans1r)) then
            rst = .false.
            print '(A)', "TEST FAILED: Beta Test 2 - 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function incomplete_gamma_test_1() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: a1 = 1.0d0
        real(real64), parameter :: x1 = 2.0d0
        real(real64), parameter :: upper1 = 0.135335283236612691894d0
        real(real64), parameter :: lower1 = 0.864664716763387308106d0
        real(real64), parameter :: a3 = 1.0d0
        real(real64), parameter :: x3 = 0.0d0
        real(real64), parameter :: upper3 = 1.0d0
        real(real64), parameter :: lower3 = 0.0d0
        real(real64), parameter :: a5 = 1.0d0
        real(real64), parameter :: x5 = 4.0d0
        real(real64), parameter :: upper5 = 0.01831563888873418029372d0
        real(real64), parameter :: lower5 = 0.9816843611112658197063d0

        ! Local Variables
        real(real64) :: dans

        ! Initialization
        rst = .true.

        ! Test 1
        dans = incomplete_gamma_upper(a1, x1)
        if (.not.is_equal(dans, upper1)) then
            rst = .false.
            print '(A)', "TEST FAILED: Incomplete Gamma 1 - 1"
        end if

        dans = incomplete_gamma_lower(a1, x1)
        if (.not.is_equal(dans, lower1)) then
            rst = .false.
            print '(A)', "TEST FAILED: Incomplete Gamma 1 - 2"
        end if

        ! Test 2
        dans = incomplete_gamma_upper(a3, x3)
        if (.not.is_equal(dans, upper3)) then
            rst = .false.
            print '(A)', "TEST FAILED: Incomplete Gamma 1 - 3"
        end if

        dans = incomplete_gamma_lower(a3, x3)
        if (.not.is_equal(dans, lower3)) then
            rst = .false.
            print '(A)', "TEST FAILED: Incomplete Gamma 1 - 4"
        end if

        ! Test 3
        dans = incomplete_gamma_upper(a5, x5)
        if (.not.is_equal(dans, upper5)) then
            rst = .false.
            print '(A)', "TEST FAILED: Incomplete Gamma 1 - 5"
        end if

        dans = incomplete_gamma_lower(a5, x5)
        if (.not.is_equal(dans, lower5)) then
            rst = .false.
            print '(A)', "TEST FAILED: Incomplete Gamma 1 - 6"
        end if
    end function

! ------------------------------------------------------------------------------
    function trimmed_mean_test_1() result(rst)
        use linalg, only : sort

        ! Arguments
        logical :: rst

        ! Parameters
        integer(int32), parameter :: n = 100
        integer(int32), parameter :: i1 = 5
        integer(int32), parameter :: i2 = n - i1 + 1
        real(real64), parameter :: p = 0.05d0

        ! Local Variables
        real(real64) :: x(n), tm, ans

        ! Initialization
        rst = .true.
        call random_number(x)
        call sort(x, .true.) ! sort into ascending order
        ans = mean(x(i1:i2))

        ! Test 1
        tm = trimmed_mean(x, p = p)
        if (.not.is_equal(tm, ans)) then
            rst = .false.
            print '(A)', "TEST FAILED: Trimmed Mean Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
end module
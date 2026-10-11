module fstats_data_column_tests
    use iso_fortran_env
    use fortran_test_helper
    use fstats_data_columns
    use strings
    implicit none
    private
    public :: test_data_columns, test_data_column_invalid_input

contains
    function test_data_columns() result(rst)
        logical :: rst, local

        rst = test_integer_column()
        local = test_real_column()
        rst = rst .and. local
        local = test_logical_column()
        rst = rst .and. local
        local = test_string_column()
        rst = rst .and. local
    end function

    subroutine check(condition, rst, label)
        logical, intent(in) :: condition
        logical, intent(inout) :: rst
        character(len = *), intent(in) :: label

        if (condition) return
        rst = .false.
        print '(A)', 'TEST FAILED: data columns - ' // label
    end subroutine

    subroutine check_values(column, expected, rst, label)
        class(data_column), intent(in) :: column
        class(*), intent(in), dimension(:) :: expected
        logical, intent(inout) :: rst
        character(len = *), intent(in) :: label
        integer(int32), allocatable, dimension(:) :: integers
        real(real64), allocatable, dimension(:) :: reals
        logical, allocatable, dimension(:) :: logicals
        type(string), allocatable, dimension(:) :: texts
        integer(int32) :: row

        call check(column%count() == size(expected), rst, label // ' count')
        if (column%count() /= size(expected)) return
        select type (expected)
        type is (integer(int32))
            allocate(integers(size(expected)))
            call column%copy(integers)
            call check(assert(integers, expected), rst, label)
        type is (real(real64))
            allocate(reals(size(expected)))
            call column%copy(reals)
            call check(assert(reals, expected, 1.0d-12), rst, label)
        type is (logical)
            allocate(logicals(size(expected)))
            call column%copy(logicals)
            call check(all(logicals .eqv. expected), rst, label)
        type is (string)
            allocate(texts(size(expected)))
            call column%copy(texts)
            do row = 1, size(expected)
                call check(char(texts(row)) == char(expected(row)), rst, label)
            end do
        class default
            call check(.false., rst, 'unsupported test type')
        end select
    end subroutine

    subroutine check_rows(column, values, rows, rst, label)
        class(data_column), intent(in) :: column
        class(*), intent(in), dimension(:) :: values
        integer(int32), intent(in), dimension(:) :: rows
        logical, intent(inout) :: rst
        character(len = *), intent(in) :: label

        select type (values)
        type is (integer(int32))
            call check_values(column, values(rows), rst, label)
        type is (real(real64))
            call check_values(column, values(rows), rst, label)
        type is (logical)
            call check_values(column, values(rows), rst, label)
        type is (string)
            call check_values(column, values(rows), rst, label)
        end select
    end subroutine

    subroutine exercise_mutations(column, values, dtype, rst)
        class(data_column), target, intent(inout) :: column
        class(*), intent(in), dimension(:) :: values
        integer(int32), intent(in) :: dtype
        logical, intent(inout) :: rst
        class(*), pointer, dimension(:) :: view

        call check(column%count() == 0, rst, 'uninitialized count')
        call check(column%get_data_type() == dtype, rst, 'data type')
        view => column%view()
        call check(.not.associated(view), rst, 'uninitialized view')
        nullify(view)
        call check_values(column, values(:0), rst, 'uninitialized copy')
        call column%set_name('original')
        call check(column%get_name() == 'original', rst, 'get name')
        call column%set_name('renamed column')
        call check(column%get_name() == 'renamed column', rst, 'rename')

        call column%append(values(2))
        call check_values(column, values(2:2), rst, 'append uninitialized')
        call column%append(values(4))
        call check_rows(column, values, [2, 4], rst, 'append populated')
        call column%insert(1, values(1))
        call check_rows(column, values, [1, 2, 4], rst, 'insert first')
        call column%insert(3, values(3))
        call check_values(column, values, rst, 'insert middle')
        call column%insert(column%count() + 1, values(1))
        call check_rows(column, values, [1, 2, 3, 4, 1], rst, 'insert end')
        call column%remove(column%count())
        call check_values(column, values, rst, 'remove last')
        call column%remove(2)
        call check_rows(column, values, [1, 3, 4], rst, 'remove middle')
        call column%remove(1)
        call check_values(column, values(3:4), rst, 'remove first')
        call column%remove(1)
        call column%remove(1)
        call check_values(column, values(:0), rst, 'remove to empty')
        view => column%view()
        call check(associated(view), rst, 'allocated empty view')
        if (associated(view)) call check(size(view) == 0, rst, 'empty view size')
        nullify(view)
        call column%insert(1, values(4))
        call check_values(column, values(4:4), rst, 'insert empty')
        call column%remove(1)
        call column%append(values(1))
        call check_values(column, values(1:1), rst, 'append empty')
        call check(column%get_name() == 'renamed column', rst, 'mutations preserve name')
        call column%set_name('')
        call check(column%get_name() == '', rst, 'empty name')
    end subroutine

    function test_integer_column() result(rst)
        logical :: rst
        type(integer_column), target :: column
        integer(int32), dimension(4) :: values
        integer(int32), dimension(2) :: copied
        class(*), pointer, dimension(:) :: view

        rst = .true.
        values = [1, 2, 3, 4]
        call exercise_mutations(column, values, FS_COLUMN_INTEGER, rst)
        column = integer_column(2)
        call check_values(column, [0, 0], rst, 'integer size constructor')
        column = integer_column(values(:2))
        values(1) = 99
        call check_values(column, [1, 2], rst, 'integer array constructor independence')
        call column%set(2, 7)
        call check(column%get(2) == 7, rst, 'integer get/set')
        call column%copy(copied)
        copied(1) = 99
        call check(column%get(1) == 1, rst, 'integer copy independence')
        view => column%view()
        select type (view)
        type is (integer(int32))
            call check(size(view) == 2, rst, 'integer view size')
            view(1) = 8
            call check(column%get(1) == 8, rst, 'integer view writes')
            call column%set(2, 9)
            call check(view(2) == 9, rst, 'integer view reads')
        class default
            call check(.false., rst, 'integer view type')
        end select
        nullify(view)
        call check(column%to_string(1) == '8', rst, 'integer default format')
        call check(column%to_string(2, 'I3.3') == '009', rst, 'integer custom format')
        column = integer_column(0)
        call check_values(column, values(:0), rst, 'integer zero constructor')
    end function

    function test_real_column() result(rst)
        logical :: rst
        type(real_column), target :: column
        real(real64), dimension(4) :: values
        real(real64), dimension(2) :: copied
        class(*), pointer, dimension(:) :: view
        character(len = :), allocatable :: text
        real(real64) :: parsed

        rst = .true.
        values = [1.0d0, 2.0d0, 3.0d0, 4.0d0]
        call exercise_mutations(column, values, FS_COLUMN_REAL, rst)
        column = real_column(2)
        call check_values(column, [0.0d0, 0.0d0], rst, 'real size constructor')
        column = real_column(values(:2))
        values(1) = 99.0d0
        call check_values(column, [1.0d0, 2.0d0], rst, 'real array constructor independence')
        call column%set(2, 7.0d0)
        call check(assert(column%get(2), 7.0d0), rst, 'real get/set')
        call column%copy(copied)
        copied(1) = 99.0d0
        call check(assert(column%get(1), 1.0d0), rst, 'real copy independence')
        view => column%view()
        select type (view)
        type is (real(real64))
            call check(size(view) == 2, rst, 'real view size')
            view(1) = 8.5d0
            call check(assert(column%get(1), 8.5d0), rst, 'real view writes')
            call column%set(2, 9.0d0)
            call check(assert(view(2), 9.0d0), rst, 'real view reads')
        class default
            call check(.false., rst, 'real view type')
        end select
        nullify(view)
        call check(column%to_string(1) == '8.5', rst, 'real default format')
        call check(column%to_string(1, 'F10.3') == '8.5', rst, 'real custom format')
        call column%set(1, 1.0d-10)
        text = column%to_string(1)
        read(text, *) parsed
        call check(assert(parsed, 1.0d-10, 1.0d-20), rst, 'real negative exponent')
        call column%set(1, 1.0d10)
        text = column%to_string(1)
        read(text, *) parsed
        call check(assert(parsed, 1.0d10, 1.0d-2), rst, 'real positive exponent')
        column = real_column(0)
        call check_values(column, values(:0), rst, 'real zero constructor')
    end function

    function test_logical_column() result(rst)
        logical :: rst
        type(logical_column), target :: column
        logical, dimension(4) :: values
        logical, dimension(2) :: copied
        class(*), pointer, dimension(:) :: view

        rst = .true.
        values = [.true., .false., .true., .false.]
        call exercise_mutations(column, values, FS_COLUMN_LOGICAL, rst)
        column = logical_column(2)
        call check_values(column, [.false., .false.], rst, 'logical size constructor')
        column = logical_column(values(:2))
        values(1) = .false.
        call check_values(column, [.true., .false.], rst, 'logical array constructor independence')
        call column%set(2, .true.)
        call check(column%get(2), rst, 'logical get/set')
        call column%copy(copied)
        copied(1) = .false.
        call check(column%get(1), rst, 'logical copy independence')
        view => column%view()
        select type (view)
        type is (logical)
            call check(size(view) == 2, rst, 'logical view size')
            view(1) = .false.
            call check(.not.column%get(1), rst, 'logical view writes')
            call column%set(2, .false.)
            call check(.not.view(2), rst, 'logical view reads')
        class default
            call check(.false., rst, 'logical view type')
        end select
        nullify(view)
        call check(column%to_string(1) == 'F', rst, 'logical default format')
        call column%set(2, .true.)
        call check(column%to_string(2, 'L2') == 'T', rst, 'logical custom format')
        column = logical_column(0)
        call check_values(column, values(:0), rst, 'logical zero constructor')
    end function

    function test_string_column() result(rst)
        logical :: rst
        type(string_column), target :: column
        type(string), dimension(4) :: values
        type(string), dimension(2) :: copied
        class(*), pointer, dimension(:) :: view

        rst = .true.
        values(1) = 'a'
        values(2) = 'different length'
        values(3) = ''
        values(4) = 'last'
        call exercise_mutations(column, values, FS_COLUMN_STRING, rst)
        column = string_column(2)
        call check(column%count() == 2, rst, 'string size constructor')
        call column%set(1, '')
        call check(column%get(1) == '', rst, 'string empty value')
        column = string_column(values(:2))
        values(1) = 'changed source'
        call check(column%get(1) == 'a', rst, 'string array constructor independence')
        call column%set(2, 'replacement')
        call check(column%get(2) == 'replacement', rst, 'string get/set')
        call column%copy(copied)
        copied(1) = 'changed copy'
        call check(column%get(1) == 'a', rst, 'string copy independence')
        view => column%view()
        select type (view)
        type is (string)
            call check(size(view) == 2, rst, 'string view size')
            view(1) = 'borrowed'
            call check(column%get(1) == 'borrowed', rst, 'string view writes')
            call column%set(2, 'visible')
            call check(char(view(2)) == 'visible', rst, 'string view reads')
        class default
            call check(.false., rst, 'string view type')
        end select
        nullify(view)
        call check(column%to_string(1) == 'borrowed', rst, 'string default format')
        call check(column%to_string(1, 'A') == 'borrowed', rst, 'string custom format')
        column = string_column(0)
        call check_values(column, values(:0), rst, 'string zero constructor')
        call column%append('character append')
        call column%insert(1, '')
        call column%insert(2, 'character insert')
        call check(column%get(1) == '', rst, 'character insert empty')
        call check(column%get(2) == 'character insert', rst, 'character insert')
        call check(column%get(3) == 'character append', rst, 'character append')
    end function

    subroutine test_data_column_invalid_input(dtype, invalid_case)
        character(len = *), intent(in) :: dtype, invalid_case
        class(data_column), allocatable :: column
        integer(int32), dimension(0) :: wrong_type
        integer(int32), dimension(2) :: integers
        real(real64), dimension(2) :: reals
        logical, dimension(2) :: logicals
        type(string), dimension(2) :: texts
        integer(int32) :: row, int_value
        real(real64) :: real_value
        logical :: logical_value
        character(len = :), allocatable :: text

        select case (dtype)
        case ('integer')
            allocate(integer_column :: column)
        case ('real')
            allocate(real_column :: column)
        case ('logical')
            allocate(logical_column :: column)
        case ('string')
            allocate(string_column :: column)
        case default
            stop 0
        end select

        if (invalid_case == 'append-type') then
            call column%append((1.0d0, 0.0d0))
        else if (invalid_case == 'insert-type') then
            call column%insert(1, (1.0d0, 0.0d0))
        else if (invalid_case == 'copy-type') then
            select case (dtype)
            case ('integer')
                call column%copy(reals(:0))
            case default
                call column%copy(wrong_type)
            end select
        else if (invalid_case == 'insert-low') then
            call column%insert(0, 1)
        else if (invalid_case == 'insert-high') then
            call column%insert(2, 1)
        else if (invalid_case == 'remove-empty') then
            call column%remove(1)
        else
            if (index(invalid_case, 'unallocated') == 0 .and. invalid_case /= 'negative-size') then
                select type (column)
                type is (integer_column)
                    call column%append(1)
                type is (real_column)
                    call column%append(1.0d0)
                type is (logical_column)
                    call column%append(.true.)
                type is (string_column)
                    call column%append('one')
                end select
            end if
            row = 1
            if (index(invalid_case, '-low') > 0) row = 0
            if (index(invalid_case, '-high') > 0) row = 2
            select case (invalid_case)
            case ('remove-low', 'remove-high')
                call column%remove(row)
            case ('format-low', 'format-high', 'format-unallocated')
                text = column%to_string(row)
            case default
                select type (column)
                type is (integer_column)
                    select case (invalid_case)
                    case ('negative-size')
                        column = integer_column(-1)
                    case ('get-low', 'get-high', 'get-unallocated')
                        int_value = column%get(row)
                    case ('set-low', 'set-high', 'set-unallocated')
                        call column%set(row, 1)
                    case ('copy-small')
                        call column%copy(integers(:0))
                    case ('copy-large')
                        call column%copy(integers)
                    end select
                type is (real_column)
                    select case (invalid_case)
                    case ('negative-size')
                        column = real_column(-1)
                    case ('get-low', 'get-high', 'get-unallocated')
                        real_value = column%get(row)
                    case ('set-low', 'set-high', 'set-unallocated')
                        call column%set(row, 1.0d0)
                    case ('copy-small')
                        call column%copy(reals(:0))
                    case ('copy-large')
                        call column%copy(reals)
                    end select
                type is (logical_column)
                    select case (invalid_case)
                    case ('negative-size')
                        column = logical_column(-1)
                    case ('get-low', 'get-high', 'get-unallocated')
                        logical_value = column%get(row)
                    case ('set-low', 'set-high', 'set-unallocated')
                        call column%set(row, .true.)
                    case ('copy-small')
                        call column%copy(logicals(:0))
                    case ('copy-large')
                        call column%copy(logicals)
                    end select
                type is (string_column)
                    select case (invalid_case)
                    case ('negative-size')
                        column = string_column(-1)
                    case ('get-low', 'get-high', 'get-unallocated')
                        text = column%get(row)
                    case ('set-low', 'set-high', 'set-unallocated')
                        call column%set(row, 'one')
                    case ('copy-small')
                        call column%copy(texts(:0))
                    case ('copy-large')
                        call column%copy(texts)
                    end select
                end select
            end select
        end if
    end subroutine
end module
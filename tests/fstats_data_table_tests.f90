module fstats_data_table_tests
    use iso_fortran_env
    use fstats_data_table
    use fstats_data_columns
    implicit none
    private
    public :: test_data_table, test_data_table_invalid_input

contains
    function test_data_table() result(rst)
        !! Exercises all table bindings, heterogeneous storage, and owned copies.
        logical :: rst
        type(data_table) :: table, empty_rows, view_source, copied_table, replacement_table
        type(data_table_view) :: live_view, single_row_view, detached_view
        type(integer_column) :: integers
        type(real_column) :: reals
        type(logical_column) :: logicals
        type(string_column) :: texts
        class(data_column), pointer :: column, second
        class(*), pointer, dimension(:) :: values
        integer(int32) :: column_index

        rst = .true.
        call check(table%column_count() == 0, rst, 'initial column count')
        call check(table%row_count() == 0, rst, 'initial row count')
        call table%remove_column(1)
        call check(table%column_count() == 0, rst, 'remove from empty table')

        integers = integer_column([1, 2, 3])
        reals = real_column([1.5d0, 2.5d0, 3.5d0])
        logicals = logical_column([.true., .false., .true.])
        texts = string_column(3)
        call texts%set(1, 'first')
        call texts%set(2, 'different length')
        call texts%set(3, '')
        call integers%set_name('integers')
        call reals%set_name('reals')
        call logicals%set_name('logicals')
        call texts%set_name('texts')

        call table%add_column(integers)
        call table%add_column(reals)
        call table%add_column(logicals)
        call table%add_column(texts)
        call check(table%column_count() == 4, rst, 'add column count')
        call check(table%row_count() == 3, rst, 'add row count')
        call integers%set(1, 99)
        call reals%set(1, 99.0d0)
        call logicals%set(1, .false.)
        call texts%set(1, 'changed source')
        call integers%set_name('changed source name')
        column => table%get_column(1)
        call check(associated(column), rst, 'integer pointer associated')
        select type (column)
        type is (integer_column)
            call check(column%get(1) == 1, rst, 'integer owned copy')
            call check(column%get_name() == 'integers', rst, 'owned name copy')
            call column%set(2, 42)
            call check(integers%get(2) == 2, rst, 'pointer writes independent of source')
        class default
            call check(.false., rst, 'integer dynamic type')
        end select
        second => table%get_column(1)
        call check(associated(column, second), rst, 'repeated access aliases stored object')
        select type (second)
        type is (integer_column)
            call check(second%get(2) == 42, rst, 'pointer writes visible in table')
        end select
        nullify(column, second)
        column => table%get_column(2)
        select type (column)
        type is (real_column)
            call check(column%get(1) == 1.5d0, rst, 'real owned copy')
        class default
            call check(.false., rst, 'real dynamic type')
        end select
        nullify(column)
        column => table%get_column(3)
        select type (column)
        type is (logical_column)
            call check(column%get(1), rst, 'logical owned copy')
        class default
            call check(.false., rst, 'logical dynamic type')
        end select
        nullify(column)
        column => table%get_column(4)
        select type (column)
        type is (string_column)
            call check(column%get(1) == 'first', rst, 'string owned copy')
            call check(column%get(2) == 'different length', rst, 'string lengths')
            call check(column%get(3) == '', rst, 'empty string')
            call column%set(1, 'table edit')
        class default
            call check(.false., rst, 'string dynamic type')
        end select
        nullify(column)
        call check(texts%get(1) == 'changed source', rst, 'string pointer write independence')

        call table%insert_column(1, reals)
        call table%insert_column(3, logicals)
        call table%insert_column(table%column_count() + 1, texts)
        call check_order(table, ['reals   ', 'integers', 'logicals', 'reals   ', &
            'logicals', 'texts   ', 'texts   '], rst, 'insert ordering')
        call check(table%row_count() == 3, rst, 'insert preserves row count')
        call reals%set(1, -1.0d0)
        column => table%get_column(1)
        select type (column)
        type is (real_column)
            call check(column%get(1) == 99.0d0, rst, 'insert owned copy')
        class default
            call check(.false., rst, 'insert dynamic type')
        end select
        nullify(column)
        call table%remove_column(3)
        call check_order(table, ['reals   ', 'integers', 'reals   ', 'logicals', &
            'texts   ', 'texts   '], rst, 'remove middle ordering')
        call table%remove_column(table%column_count())
        call check_order(table, ['reals   ', 'integers', 'reals   ', 'logicals', &
            'texts   '], rst, 'remove last ordering')
        call table%remove_column(1)
        call check_order(table, ['integers', 'reals   ', 'logicals', 'texts   '], rst, 'remove first ordering')
        do column_index = 1, 4
            call table%remove_column(1)
        end do
        call check(table%column_count() == 0, rst, 'remove all column count')
        call check(table%row_count() == 0, rst, 'remove all row count')

        call table%insert_column(1, integer_column([8, 9]))
        call check(table%row_count() == 2, rst, 'insert into empty sets row count')
        call table%remove_column(1)
        call table%add_column(real_column([7.0d0]))
        call check(table%row_count() == 1, rst, 'reuse table with different row count')
        call add_local_column(table)
        column => table%get_column(2)
        select type (column)
        type is (string_column)
            call check(column%get(1) == 'local value', rst, 'source lifetime independence')
        class default
            call check(.false., rst, 'local column dynamic type')
        end select
        nullify(column)

        call empty_rows%add_column(integer_column(0))
        call empty_rows%insert_column(1, string_column(0))
        call check(empty_rows%column_count() == 2, rst, 'zero-row columns stored')
        call check(empty_rows%row_count() == 0, rst, 'zero-row table count')
        column => empty_rows%get_column(1)
        call check(column%count() == 0, rst, 'get zero-row column')
        nullify(column)

        call view_source%add_column(integer_column([10, 20, 30, 40]))
        call view_source%add_column(real_column([1.0d0, 2.0d0, 3.0d0, 4.0d0]))
        live_view = view_source%view(2, 3, 2, 2)
        call check(live_view%row_count() == 2, rst, 'view initial row count')
        call check(live_view%column_count() == 1, rst, 'view initial column count')
        values => live_view%column_slice(1)
        select type (values)
        type is (real(real64))
            call check(values(1) == 2.0d0 .and. values(2) == 3.0d0, &
                rst, 'view selects requested rows')
            values(1) = 25.0d0
        class default
            call check(.false., rst, 'view real slice dynamic type')
        end select
        nullify(values)
        column => view_source%get_column(2)
        select type (column)
        type is (real_column)
            call check(column%get(2) == 25.0d0, rst, 'view writes shared source data')
        class default
            call check(.false., rst, 'view real source dynamic type')
        end select
        nullify(column)

        call view_source%insert_column(1, integer_column([5, 6, 7, 8]))
        values => live_view%column_slice(1)
        select type (values)
        type is (integer(int32))
            call check(values(1) == 20 .and. values(2) == 30, &
                rst, 'view follows current positional columns after insertion')
        class default
            call check(.false., rst, 'view inserted integer slice dynamic type')
        end select
        nullify(values)

        single_row_view = view_source%slice(2, 2)
        call check(single_row_view%row_count() == 1, rst, 'single-row slice count')
        values => single_row_view%column_slice(1)
        select type (values)
        type is (integer(int32))
            call check(size(values) == 1 .and. values(1) == 6, &
                rst, 'single-row slice value')
        class default
            call check(.false., rst, 'single-row slice dynamic type')
        end select
        nullify(values)

        call view_source%remove_column(2)
        values => live_view%column_slice(1)
        select type (values)
        type is (real(real64))
            call check(values(1) == 25.0d0 .and. values(2) == 3.0d0, &
                rst, 'view follows current positional columns after removal')
        class default
            call check(.false., rst, 'view removed real slice dynamic type')
        end select
        nullify(values)
        call view_source%remove_column(2)
        call check(live_view%column_count() == 0, rst, 'view remains valid when selection disappears')
        call view_source%add_column(real_column([9.0d0, 10.0d0, 11.0d0, 12.0d0]))
        call check(live_view%column_count() == 1, rst, 'view reflects later positional insertion')

        call replacement_table%add_column(integer_column([100, 200, 300, 400]))
        call replacement_table%add_column(real_column([110.0d0, 210.0d0, 310.0d0, 410.0d0]))
        view_source = replacement_table
        values => live_view%column_slice(1)
        select type (values)
        type is (real(real64))
            call check(values(1) == 210.0d0 .and. values(2) == 310.0d0, &
                rst, 'view remains attached across table assignment')
        class default
            call check(.false., rst, 'assigned view slice dynamic type')
        end select
        nullify(values)

        detached_view = make_detached_view()
        call check(detached_view%row_count() == 2, rst, 'view retains state after table destruction')
        values => detached_view%column_slice(1)
        select type (values)
        type is (integer(int32))
            call check(values(1) == 12 .and. values(2) == 13, &
                rst, 'detached view retains values')
        class default
            call check(.false., rst, 'detached view slice dynamic type')
        end select
        nullify(values)

        copied_table = view_source
        column => copied_table%get_column(1)
        select type (column)
        type is (integer_column)
            call column%set(1, -1)
        class default
            call check(.false., rst, 'copied table dynamic type')
        end select
        nullify(column)
        column => view_source%get_column(1)
        select type (column)
        type is (integer_column)
            call check(column%get(1) == 100, rst, 'original table dynamic type after copy')
        class default
            call check(.false., rst, 'original table dynamic type after copy')
        end select
        nullify(column)
    end function

    function make_detached_view() result(rst)
        !! Creates a view whose source table is finalized before the caller uses it.
        type(data_table_view) :: rst
            !! The view retaining the local table's shared state.

        type(data_table) :: source

        call source%add_column(integer_column([11, 12, 13]))
        rst = source%view(2, 3, 1, 1)
    end function

    subroutine add_local_column(table)
        !! Ensures the table owns its copy after a source column goes out of scope.
        class(data_table), intent(inout) :: table
        type(string_column) :: local

        local = string_column(1)
        call local%set(1, 'local value')
        call table%add_column(local)
    end subroutine

    subroutine check_order(table, names, rst, label)
        class(data_table), intent(in) :: table
        character(len = *), intent(in), dimension(:) :: names
        logical, intent(inout) :: rst
        character(len = *), intent(in) :: label
        class(data_column), pointer :: column
        integer(int32) :: column_index

        call check(table%column_count() == size(names), rst, label // ' count')
        if (table%column_count() /= size(names)) return
        do column_index = 1, size(names)
            column => table%get_column(column_index)
            call check(column%get_name() == trim(names(column_index)), rst, label)
        end do
        nullify(column)
    end subroutine

    subroutine check(condition, rst, label)
        logical, intent(in) :: condition
        logical, intent(inout) :: rst
        character(len = *), intent(in) :: label

        if (condition) return
        rst = .false.
        print '(A)', 'TEST FAILED: data table - ' // label
    end subroutine

    subroutine test_data_table_invalid_input(invalid_case)
        !! Exercises bounds and row-count validation in separate processes.
        character(len = *), intent(in) :: invalid_case
        type(data_table) :: table
        class(data_column), pointer :: column

        if (index(invalid_case, 'empty') == 0) call table%add_column(integer_column([1, 2]))
        select case (invalid_case)
        case ('get-empty', 'get-high')
            column => table%get_column(table%column_count() + 1)
        case ('get-low')
            column => table%get_column(0)
        case ('remove-low')
            call table%remove_column(0)
        case ('remove-high')
            call table%remove_column(table%column_count() + 1)
        case ('insert-low', 'insert-empty-low')
            call table%insert_column(0, integer_column([1, 2]))
        case ('insert-high', 'insert-empty-high')
            call table%insert_column(table%column_count() + 2, integer_column([1, 2]))
        case ('add-short')
            call table%add_column(real_column(1))
        case ('add-long')
            call table%add_column(real_column(3))
        case ('insert-short')
            call table%insert_column(1, real_column(1))
        case ('insert-long')
            call table%insert_column(1, real_column(3))
        case ('add-zero')
            call table%add_column(real_column(0))
        case ('insert-zero')
            call table%insert_column(1, real_column(0))
        case ('zero-add-nonempty')
            call table%remove_column(1)
            call table%add_column(integer_column(0))
            call table%add_column(real_column(1))
        case ('zero-insert-nonempty')
            call table%remove_column(1)
            call table%add_column(integer_column(0))
            call table%insert_column(1, real_column(1))
        end select
    end subroutine
end module
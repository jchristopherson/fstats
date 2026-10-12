module fstats_data_table_io_tests
    use iso_fortran_env
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use fstats_data_table
    use fstats_data_columns
    use fstats_data_table_io
    implicit none
    private
    public :: test_data_table_io, test_data_table_io_invalid_input

    character, parameter :: LF = achar(10)
    character, parameter :: CR = achar(13)
    character, parameter :: TAB = achar(9)

contains
    function test_data_table_io() result(rst)
        !! Exercises delimited file reading, type inference, and round trips.
        logical :: rst

        rst = .true.
        call test_inference(rst)
        call test_header_row(rst)
        call test_no_header(rst)
        call test_quoting_and_edge_cases(rst)
        call test_round_trip(rst)
    end function

! ------------------------------------------------------------------------------
    subroutine test_inference(rst)
        !! Default comma delimiter, first-row header, CRLF endings, blank rows,
        !! missing values, quoting, and inference of every column type.
        logical, intent(inout) :: rst
        character(len = *), parameter :: fname = "fstats_io_inference.csv"
        type(data_table) :: table
        class(data_column), pointer :: column

        call write_text(fname, &
            "id,value,count,flag,name,empty,partial_flag,mixed" // CR // LF // &
            "1,1.5,10,TRUE,""Smith, John"",,T,1" // CR // LF // &
            "2, 2 ,,false,""He said """"hi"""""",,,T" // CR // LF // &
            CR // LF // &
            "   " // TAB // CR // LF // &
            "3,-3.25e2,30,.T.,  plain text  ,,F,2" // CR // LF)
        call read_delimited_table(fname, table)
        call delete_file(fname)

        call check(table%column_count() == 8, rst, 'inference column count')
        call check(table%row_count() == 3, rst, 'inference skips blank rows')
        if (table%column_count() /= 8 .or. table%row_count() /= 3) return

        column => table%get_column(1)
        call check(column%get_name() == 'id', rst, 'inference header name')
        select type (column)
        type is (integer_column)
            call check(all([column%get(1), column%get(2), column%get(3)] == [1, 2, 3]), &
                rst, 'integer values')
        class default
            call check(.false., rst, 'integer column inferred')
        end select

        column => table%get_column(2)
        select type (column)
        type is (real_column)
            call check(column%get(1) == 1.5d0 .and. column%get(2) == 2.0d0 .and. &
                column%get(3) == -325.0d0, rst, 'mixed integer and real values')
        class default
            call check(.false., rst, 'mixed integer and real inferred as real')
        end select

        column => table%get_column(3)
        select type (column)
        type is (real_column)
            call check(column%get(1) == 10.0d0 .and. ieee_is_nan(column%get(2)) .and. &
                column%get(3) == 30.0d0, rst, 'missing integer stored as NaN')
        class default
            call check(.false., rst, 'integer with missing inferred as real')
        end select

        column => table%get_column(4)
        select type (column)
        type is (logical_column)
            call check(column%get(1) .and. .not.column%get(2) .and. column%get(3), &
                rst, 'logical values')
        class default
            call check(.false., rst, 'logical column inferred')
        end select

        column => table%get_column(5)
        select type (column)
        type is (string_column)
            call check(column%get(1) == 'Smith, John', rst, 'quoted delimiter')
            call check(column%get(2) == 'He said "hi"', rst, 'escaped quotes')
            call check(column%get(3) == 'plain text', rst, 'unquoted text trimmed')
        class default
            call check(.false., rst, 'string column inferred')
        end select

        column => table%get_column(6)
        select type (column)
        type is (real_column)
            call check(all(ieee_is_nan([column%get(1), column%get(2), column%get(3)])), &
                rst, 'empty column is NaN')
        class default
            call check(.false., rst, 'empty column inferred as real')
        end select

        column => table%get_column(7)
        select type (column)
        type is (string_column)
            call check(column%get(1) == 'T' .and. column%get(2) == '' .and. &
                column%get(3) == 'F', rst, 'logical with missing kept as text')
        class default
            call check(.false., rst, 'logical with missing inferred as string')
        end select

        column => table%get_column(8)
        select type (column)
        type is (string_column)
            call check(column%get(2) == 'T', rst, 'mixed logical and numeric text')
        class default
            call check(.false., rst, 'mixed logical and numeric inferred as string')
        end select
    end subroutine

! ------------------------------------------------------------------------------
    subroutine test_header_row(rst)
        !! A header below a preamble, with a semicolon delimiter.
        logical, intent(inout) :: rst
        character(len = *), parameter :: fname = "fstats_io_header_row.csv"
        type(data_table) :: table
        class(data_column), pointer :: column

        call write_text(fname, &
            "# instrument log" // LF // &
            "units;s;V;extra" // LF // &
            "time;volts" // LF // &
            "0.0;1" // LF // &
            "0.5;2" // LF)
        call read_delimited_table(fname, table, delimiter = ';', header_row = 3)
        call delete_file(fname)

        call check(table%column_count() == 2, rst, 'header row column count')
        call check(table%row_count() == 2, rst, 'header row row count')
        if (table%column_count() /= 2) return
        column => table%get_column(1)
        call check(column%get_name() == 'time', rst, 'header row first name')
        select type (column)
        type is (real_column)
            call check(column%get(2) == 0.5d0, rst, 'header row real value')
        class default
            call check(.false., rst, 'header row real column')
        end select
        column => table%get_column(2)
        call check(column%get_name() == 'volts', rst, 'header row second name')
        select type (column)
        type is (integer_column)
            call check(column%get(2) == 2, rst, 'header row integer value')
        class default
            call check(.false., rst, 'header row integer column')
        end select
    end subroutine

! ------------------------------------------------------------------------------
    subroutine test_no_header(rst)
        !! No header, tab delimiter, ragged rows, and no final line terminator.
        logical, intent(inout) :: rst
        character(len = *), parameter :: fname = "fstats_io_no_header.txt"
        type(data_table) :: table
        class(data_column), pointer :: column

        call write_text(fname, &
            "1" // TAB // "2.5" // TAB // "x" // LF // &
            "3" // LF // LF // &
            "4" // TAB // TAB // "y")
        call read_delimited_table(fname, table, delimiter = TAB, header_row = 0)
        call delete_file(fname)

        call check(table%column_count() == 3, rst, 'no header column count')
        call check(table%row_count() == 3, rst, 'no header row count')
        if (table%column_count() /= 3 .or. table%row_count() /= 3) return
        column => table%get_column(1)
        call check(column%get_name() == 'column_1', rst, 'default column name')
        select type (column)
        type is (integer_column)
            call check(column%get(3) == 4, rst, 'no header integer value')
        class default
            call check(.false., rst, 'no header integer column')
        end select
        column => table%get_column(2)
        call check(column%get_name() == 'column_2', rst, 'second default column name')
        select type (column)
        type is (real_column)
            call check(column%get(1) == 2.5d0 .and. ieee_is_nan(column%get(2)) .and. &
                ieee_is_nan(column%get(3)), rst, 'ragged and empty fields are NaN')
        class default
            call check(.false., rst, 'no header real column')
        end select
        column => table%get_column(3)
        select type (column)
        type is (string_column)
            call check(column%get(1) == 'x' .and. column%get(2) == '' .and. &
                column%get(3) == 'y', rst, 'ragged string padded')
        class default
            call check(.false., rst, 'no header string column')
        end select
    end subroutine

! ------------------------------------------------------------------------------
    subroutine test_quoting_and_edge_cases(rst)
        !! Byte order marks, quoted line breaks, header-only and empty files,
        !! large integers, and special real values.
        logical, intent(inout) :: rst
        character(len = *), parameter :: fname = "fstats_io_edge.csv"
        type(data_table) :: table
        class(data_column), pointer :: column

        call write_text(fname, achar(239) // achar(187) // achar(191) // &
            "a,b,big,special" // LF // &
            """line 1" // LF // "line 2"",3,3000000000,NaN" // LF // &
            """"",4,1,-Inf" // LF)
        call read_delimited_table(fname, table)
        call check(table%column_count() == 4, rst, 'edge column count')
        call check(table%row_count() == 2, rst, 'quoted line break row count')
        if (table%column_count() == 4 .and. table%row_count() == 2) then
            column => table%get_column(1)
            call check(column%get_name() == 'a', rst, 'byte order mark skipped')
            select type (column)
            type is (string_column)
                call check(column%get(1) == 'line 1' // LF // 'line 2', rst, &
                    'quoted line break')
                call check(column%get(2) == '', rst, 'quoted empty string')
            class default
                call check(.false., rst, 'quoted line break column type')
            end select
            column => table%get_column(3)
            select type (column)
            type is (real_column)
                call check(column%get(1) == 3.0d9, rst, 'int32 overflow read as real')
            class default
                call check(.false., rst, 'int32 overflow column type')
            end select
            column => table%get_column(4)
            select type (column)
            type is (real_column)
                call check(ieee_is_nan(column%get(1)) .and. column%get(2) < -huge(1.0d0), &
                    rst, 'NaN and infinity')
            class default
                call check(.false., rst, 'special values column type')
            end select
        end if

        call write_text(fname, "x,y" // LF)
        call read_delimited_table(fname, table)
        call check(table%column_count() == 2, rst, 'header only column count')
        call check(table%row_count() == 0, rst, 'header only row count')

        call write_text(fname, "")
        call read_delimited_table(fname, table, header_row = 0)
        call check(table%column_count() == 0, rst, 'empty file column count')
        call delete_file(fname)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine test_round_trip(rst)
        !! Writes every column type with a custom delimiter and reads it back.
        logical, intent(inout) :: rst
        character(len = *), parameter :: fname = "fstats_io_round_trip.txt"
        type(data_table) :: table, copy
        type(integer_column) :: integers
        type(real_column) :: reals
        type(logical_column) :: logicals
        type(string_column) :: texts
        class(data_column), pointer :: column
        real(real64) :: nan
        real(real64), parameter :: values(3) = [0.1d0, 0.0d0, 12345.678901234567d0]

        nan = ieee_value(nan, ieee_quiet_nan)
        integers = integer_column([1, -2, 3])
        reals = real_column([values(1), nan, values(3)])
        logicals = logical_column([.true., .false., .true.])
        texts = string_column(3)
        call texts%set(1, 'a|b')
        call texts%set(2, ' padded')
        call texts%set(3, 'quote"d')
        call integers%set_name('int')
        call reals%set_name('real|value')
        call logicals%set_name('flag')
        call texts%set_name('text')
        call table%add_column(integers)
        call table%add_column(reals)
        call table%add_column(logicals)
        call table%add_column(texts)

        call write_delimited_table(fname, table, delimiter = '|', real_fmt = 'ES26.17E3')
        call read_delimited_table(fname, copy, delimiter = '|')
        call check(copy%column_count() == 4, rst, 'round trip column count')
        call check(copy%row_count() == 3, rst, 'round trip row count')
        if (copy%column_count() == 4 .and. copy%row_count() == 3) then
            column => copy%get_column(1)
            call check(column%get_name() == 'int', rst, 'round trip name')
            select type (column)
            type is (integer_column)
                call check(column%get(2) == -2, rst, 'round trip integer')
            class default
                call check(.false., rst, 'round trip integer type')
            end select
            column => copy%get_column(2)
            call check(column%get_name() == 'real|value', rst, 'round trip quoted name')
            select type (column)
            type is (real_column)
                call check(column%get(1) == values(1) .and. column%get(3) == values(3), &
                    rst, 'round trip full precision reals')
                call check(ieee_is_nan(column%get(2)), rst, 'round trip NaN')
            class default
                call check(.false., rst, 'round trip real type')
            end select
            column => copy%get_column(3)
            select type (column)
            type is (logical_column)
                call check(column%get(1) .and. .not.column%get(2), rst, 'round trip logical')
            class default
                call check(.false., rst, 'round trip logical type')
            end select
            column => copy%get_column(4)
            select type (column)
            type is (string_column)
                call check(column%get(1) == 'a|b' .and. column%get(2) == ' padded' .and. &
                    column%get(3) == 'quote"d', rst, 'round trip strings')
            class default
                call check(.false., rst, 'round trip string type')
            end select
        end if

        call write_delimited_table(fname, table, write_header = .false.)
        call read_delimited_table(fname, copy, header_row = 0)
        call delete_file(fname)
        call check(copy%column_count() == 4, rst, 'headerless round trip column count')
        call check(copy%row_count() == 3, rst, 'headerless round trip row count')
        if (copy%column_count() == 4) then
            column => copy%get_column(2)
            select type (column)
            type is (real_column)
                call check(abs(column%get(3) - values(3)) < 1.0d-6 * values(3), rst, &
                    'default real format')
            class default
                call check(.false., rst, 'headerless real type')
            end select
        end if
    end subroutine

! ------------------------------------------------------------------------------
    subroutine write_text(fname, text)
        !! Writes text verbatim so tests control line terminators exactly.
        character(len = *), intent(in) :: fname, text
        integer :: unit

        open(newunit = unit, file = fname, access = 'stream', form = 'unformatted', &
            status = 'replace', action = 'write')
        if (len(text) > 0) write(unit) text
        close(unit)
    end subroutine

    subroutine delete_file(fname)
        character(len = *), intent(in) :: fname
        integer :: unit, ios

        open(newunit = unit, file = fname, status = 'old', iostat = ios)
        if (ios == 0) close(unit, status = 'delete')
    end subroutine

    subroutine check(condition, rst, label)
        logical, intent(in) :: condition
        logical, intent(inout) :: rst
        character(len = *), intent(in) :: label

        if (condition) return
        rst = .false.
        print '(A)', 'TEST FAILED: data table io - ' // label
    end subroutine

! ------------------------------------------------------------------------------
    subroutine test_data_table_io_invalid_input(invalid_case)
        !! Exercises error stops in separate processes.
        character(len = *), intent(in) :: invalid_case
        character(len = *), parameter :: fname = 'fstats_io_invalid_' // 'case.csv'
        type(data_table) :: table

        select case (invalid_case)
        case ('missing-file')
            call read_delimited_table('fstats_io_file_that_does_not_exist.csv', table)
        case ('quote-delimiter')
            call write_text(fname, 'a' // LF // '1' // LF)
            call read_delimited_table(fname, table, delimiter = '"')
        case ('negative-header')
            call write_text(fname, 'a' // LF // '1' // LF)
            call read_delimited_table(fname, table, header_row = -1)
        case ('header-past-end')
            call write_text(fname, 'a' // LF // '1' // LF)
            call read_delimited_table(fname, table, header_row = 3)
        case ('blank-header')
            call write_text(fname, LF // 'a' // LF // '1' // LF)
            call read_delimited_table(fname, table)
        case ('empty-file')
            call write_text(fname, '')
            call read_delimited_table(fname, table)
        case ('extra-fields')
            call write_text(fname, 'a,b' // LF // '1,2,3' // LF)
            call read_delimited_table(fname, table)
        case ('unterminated-quote')
            call write_text(fname, 'a,b' // LF // '1,"open' // LF)
            call read_delimited_table(fname, table)
        case ('malformed-quote')
            call write_text(fname, 'a,b' // LF // '1,"x"y' // LF)
            call read_delimited_table(fname, table)
        case ('write-quote-delimiter')
            call table%add_column(integer_column([1, 2]))
            call write_delimited_table(fname, table, delimiter = '"')
        case ('write-bad-path')
            call table%add_column(integer_column([1, 2]))
            call write_delimited_table('fstats_missing_directory/out.csv', table)
        end select
    end subroutine
end module

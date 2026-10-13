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

module fstats_data_table_io
    !! Provides routines for reading and writing [[data_table]] objects from
    !! and to delimited text files (e.g. CSV files).
    !!
    !! Both routines touch the disk exactly once: the reader loads the entire
    !! file with a single stream read and parses it in memory, and the writer
    !! assembles the entire file in memory before a single stream write.
    use iso_fortran_env
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, &
        ieee_positive_inf, ieee_negative_inf
    use strings, only : string, assignment(=)
    use fstats_errors
    use fstats_data_columns
    use fstats_data_table
    implicit none
    private
    public :: read_delimited_table
    public :: write_delimited_table

    character, parameter :: LF = achar(10)
    character, parameter :: CR = achar(13)
    character, parameter :: TAB = achar(9)
    character, parameter :: QUOTE = '"'

    integer(int32), parameter :: FIELD_LOGICAL = 1
        !! The field holds a logical value.
    integer(int32), parameter :: FIELD_INTEGER = 2
        !! The field holds a 32-bit integer value.
    integer(int32), parameter :: FIELD_REAL = 3
        !! The field holds a floating-point value.
    integer(int32), parameter :: FIELD_TEXT = 4
        !! The field holds a value that is not logical or numeric.

    type delimited_fields
        !! Describes the fields parsed from an in-memory file buffer.  Fields
        !! are stored as character extents into the buffer to avoid copying.
        integer(int32) :: field_count = 0
            !! The number of stored fields.
        integer(int32) :: row_count = 0
            !! The number of stored data rows.
        integer(int32), allocatable, dimension(:) :: first
            !! The buffer index of the first character of each field.
        integer(int32), allocatable, dimension(:) :: last
            !! The buffer index of the last character of each field; a value
            !! less than the matching first index denotes an empty field.
        logical, allocatable, dimension(:) :: quoted
            !! True for fields that were enclosed in double quotes.
        integer(int32), allocatable, dimension(:) :: row_first
            !! The index of the first field of each data row.
        integer(int32), allocatable, dimension(:) :: row_size
            !! The number of fields in each data row.
    end type

    type column_pointer
        !! Wraps a borrowed column pointer so columns can be held in an array.
        class(data_column), pointer :: ptr => null()
    end type

contains
! ******************************************************************************
! READING
! ------------------------------------------------------------------------------
    subroutine read_delimited_table(fname, table, delimiter, header_row)
        !! Reads a delimited text file into a [[data_table]].
        !!
        !! The file is read from disk with a single operation and parsed in
        !! memory.  Rows may end with LF, CRLF, or CR, and a leading UTF-8 byte
        !! order mark is ignored.  Rows containing only spaces and tabs are
        !! skipped.  Unquoted fields are trimmed of surrounding spaces and tabs
        !! (unless the delimiter is itself a space or tab).  Fields may be
        !! enclosed in double quotes to embed delimiters, line breaks, or
        !! quote characters (written as two consecutive quotes); quoting does
        !! not affect type inference.  Consecutive delimiters denote an empty,
        !! or missing, field.  Rows with fewer fields than the table has
        !! columns are padded with missing values.
        !!
        !! The type of each column is inferred from its non-empty entries:
        !!
        !! - [[logical_column]]: every entry is one of T, F, TRUE, FALSE, .T.,
        !!   .F., .TRUE., or .FALSE. (case-insensitive) and nothing is missing.
        !!
        !! - [[integer_column]]: every entry is an integer representable as a
        !!   32-bit integer and nothing is missing.
        !!
        !! - [[real_column]]: every entry is numeric (including NaN, Inf, and
        !!   Infinity), integers are mixed with reals or are missing values,
        !!   or the column is entirely empty.  Missing entries are stored as
        !!   NaN.
        !!
        !! - [[string_column]]: any other column, including logical columns
        !!   with missing entries and columns mixing logical and numeric
        !!   entries.  Missing entries are stored as empty strings.
        !!
        !! The following errors terminate the program via error stop:
        !!
        !! - [[FS_INVALID_INPUT_ERROR]]: the file cannot be opened or read,
        !!   the delimiter is invalid, header_row is negative, the header row
        !!   is blank or beyond the end of the file, or a quoted field is
        !!   unterminated or followed by text other than a delimiter.
        !!
        !! - [[FS_ARRAY_SIZE_ERROR]]: a data row holds more fields than the
        !!   header row.
        !!
        !! - [[FS_MEMORY_ERROR]]: the file exceeds 2^31 - 1 bytes.
        character(len = *), intent(in) :: fname
            !! The name of the file to read.
        type(data_table), intent(out) :: table
            !! The resulting [[data_table]], holding one column per field.
        character, intent(in), optional :: delimiter
            !! The field delimiter.  The default is a comma.  The delimiter
            !! must not be a double quote, carriage return, or line feed.
        integer(int32), intent(in), optional :: header_row
            !! The one-based row, counting blank rows, holding the column
            !! names.  Rows above the header are ignored, and data begins on
            !! the row after the header.  The default is 1.  Use 0 for files
            !! without a header; the columns are then named column_1,
            !! column_2, and so on, and the column count is the largest number
            !! of fields found in any row.

        character :: delim
        character(len = :), allocatable :: buffer, name
        character(len = 16) :: index_text
        integer(int32) :: hdr, unit, ios, start, ncols, header_first, &
            header_size, j
        integer(int64) :: fsize
        type(delimited_fields) :: fields

        delim = ","
        if (present(delimiter)) delim = delimiter
        if (delim == QUOTE .or. delim == CR .or. delim == LF) &
            error stop FS_INVALID_INPUT_ERROR
        hdr = 1
        if (present(header_row)) hdr = header_row
        if (hdr < 0) error stop FS_INVALID_INPUT_ERROR

        open(newunit = unit, file = fname, access = "stream", &
            form = "unformatted", action = "read", status = "old", &
            iostat = ios)
        if (ios /= 0) error stop FS_INVALID_INPUT_ERROR
        inquire(unit = unit, size = fsize)
        if (fsize < 0 .or. fsize > huge(1_int32)) then
            close(unit)
            error stop FS_MEMORY_ERROR
        end if
        allocate(character(len = int(fsize, int32)) :: buffer)
        if (fsize > 0) read(unit, iostat = ios) buffer
        close(unit)
        if (ios /= 0) error stop FS_INVALID_INPUT_ERROR

        start = 1
        if (len(buffer) >= 3) then
            if (buffer(1:3) == achar(239) // achar(187) // achar(191)) start = 4
        end if

        call parse_buffer(buffer, start, delim, hdr, fields, header_first, &
            header_size)

        if (hdr > 0) then
            ncols = header_size
            if (fields%row_count > 0) then
                if (maxval(fields%row_size(1:fields%row_count)) > ncols) &
                    error stop FS_ARRAY_SIZE_ERROR
            end if
        else if (fields%row_count > 0) then
            ncols = maxval(fields%row_size(1:fields%row_count))
        else
            ncols = 0
        end if

        do j = 1, ncols
            if (hdr > 0) then
                name = get_field_text(buffer, fields, header_first + j - 1)
            else
                write(index_text, "(I0)") j
                name = "column_" // trim(index_text)
            end if
            call add_inferred_column(table, buffer, fields, j, name)
        end do
    end subroutine

! ------------------------------------------------------------------------------
    subroutine parse_buffer(buffer, start, delim, header_row, fields, &
        header_first, header_size)
        !! Splits the in-memory file contents into fields in a single pass.
        character(len = *), intent(in) :: buffer
            !! The file contents.
        integer(int32), intent(in) :: start
            !! The buffer index at which parsing begins.
        character, intent(in) :: delim
            !! The field delimiter.
        integer(int32), intent(in) :: header_row
            !! The header row index, or 0 if there is no header.
        type(delimited_fields), intent(out) :: fields
            !! The parsed header and data fields.
        integer(int32), intent(out) :: header_first
            !! The index of the first header field.
        integer(int32), intent(out) :: header_size
            !! The number of header fields; zero if no header was found.

        integer(int32) :: n, pos, k, rec, first
        logical :: blank

        n = len(buffer)
        pos = start
        rec = 0
        header_first = 0
        header_size = 0
        allocate(fields%first(1024), fields%last(1024), fields%quoted(1024))
        allocate(fields%row_first(256), fields%row_size(256))

        do while (pos <= n)
            rec = rec + 1
            k = pos
            do while (k <= n)
                if (buffer(k:k) /= " " .and. buffer(k:k) /= TAB) exit
                k = k + 1
            end do
            blank = k > n
            if (.not.blank) blank = buffer(k:k) == LF .or. buffer(k:k) == CR
            if (blank) then
                if (rec == header_row) error stop FS_INVALID_INPUT_ERROR
                pos = k
                call skip_terminator(buffer, pos)
                cycle
            end if

            first = fields%field_count + 1
            call parse_record(buffer, delim, pos, fields)
            if (rec < header_row) then
                ! Discard rows above the header; they are parsed only so that
                ! quoted line breaks are tracked correctly.
                fields%field_count = first - 1
            else if (rec == header_row) then
                header_first = first
                header_size = fields%field_count - first + 1
            else
                call push_row(fields, first, fields%field_count - first + 1)
            end if
        end do

        if (header_row > 0 .and. header_size == 0) &
            error stop FS_INVALID_INPUT_ERROR
    end subroutine

! ------------------------------------------------------------------------------
    subroutine parse_record(buffer, delim, pos, fields)
        !! Parses one row of fields, advancing pos past the row terminator.
        character(len = *), intent(in) :: buffer
            !! The file contents.
        character, intent(in) :: delim
            !! The field delimiter.
        integer(int32), intent(inout) :: pos
            !! On input, the start of the row; on output, the start of the
            !! next row.
        type(delimited_fields), intent(inout) :: fields
            !! The field collection to which the row's fields are appended.

        integer(int32) :: n, s, e
        logical :: quoted
        character :: c

        n = len(buffer)
        do
            do while (pos <= n)
                if (.not.is_padding(buffer(pos:pos), delim)) exit
                pos = pos + 1
            end do

            quoted = .false.
            if (pos <= n) quoted = buffer(pos:pos) == QUOTE
            if (quoted) then
                pos = pos + 1
                s = pos
                do
                    if (pos > n) error stop FS_INVALID_INPUT_ERROR
                    if (buffer(pos:pos) == QUOTE) then
                        if (pos == n) exit
                        if (buffer(pos + 1:pos + 1) /= QUOTE) exit
                        pos = pos + 1
                    end if
                    pos = pos + 1
                end do
                e = pos - 1
                pos = pos + 1
                do while (pos <= n)
                    if (.not.is_padding(buffer(pos:pos), delim)) exit
                    pos = pos + 1
                end do
                if (pos <= n) then
                    c = buffer(pos:pos)
                    if (c /= delim .and. c /= LF .and. c /= CR) &
                        error stop FS_INVALID_INPUT_ERROR
                end if
            else
                s = pos
                do while (pos <= n)
                    c = buffer(pos:pos)
                    if (c == delim .or. c == LF .or. c == CR) exit
                    pos = pos + 1
                end do
                e = pos - 1
                do while (e >= s)
                    if (.not.is_padding(buffer(e:e), delim)) exit
                    e = e - 1
                end do
            end if

            call push_field(fields, s, e, quoted)
            if (pos > n) return
            if (buffer(pos:pos) /= delim) exit
            pos = pos + 1
        end do
        call skip_terminator(buffer, pos)
    end subroutine

! ------------------------------------------------------------------------------
    pure subroutine skip_terminator(buffer, pos)
        !! Advances past an LF, CRLF, or CR row terminator located at pos.
        character(len = *), intent(in) :: buffer
            !! The file contents.
        integer(int32), intent(inout) :: pos
            !! The buffer position.

        if (pos > len(buffer)) return
        if (buffer(pos:pos) == CR) then
            pos = pos + 1
            if (pos <= len(buffer)) then
                if (buffer(pos:pos) == LF) pos = pos + 1
            end if
        else if (buffer(pos:pos) == LF) then
            pos = pos + 1
        end if
    end subroutine

! ------------------------------------------------------------------------------
    pure function is_padding(c, delim) result(rst)
        !! Determines if a character is insignificant whitespace.
        character, intent(in) :: c
            !! The character to test.
        character, intent(in) :: delim
            !! The field delimiter, which is never treated as whitespace.
        logical :: rst
            !! True if the character is a space or tab other than the delimiter.

        rst = (c == " " .or. c == TAB) .and. c /= delim
    end function

! ------------------------------------------------------------------------------
    pure subroutine push_field(fields, first, last, quoted)
        !! Appends a field extent, growing storage geometrically as needed.
        type(delimited_fields), intent(inout) :: fields
            !! The field collection.
        integer(int32), intent(in) :: first
            !! The buffer index of the field's first character.
        integer(int32), intent(in) :: last
            !! The buffer index of the field's last character.
        logical, intent(in) :: quoted
            !! True if the field was quoted.

        integer(int32), allocatable, dimension(:) :: itemp
        logical, allocatable, dimension(:) :: ltemp
        integer(int32) :: n

        n = fields%field_count
        if (n == size(fields%first)) then
            allocate(itemp(2 * n))
            itemp(1:n) = fields%first
            call move_alloc(itemp, fields%first)
            allocate(itemp(2 * n))
            itemp(1:n) = fields%last
            call move_alloc(itemp, fields%last)
            allocate(ltemp(2 * n))
            ltemp(1:n) = fields%quoted
            call move_alloc(ltemp, fields%quoted)
        end if
        n = n + 1
        fields%first(n) = first
        fields%last(n) = last
        fields%quoted(n) = quoted
        fields%field_count = n
    end subroutine

! ------------------------------------------------------------------------------
    pure subroutine push_row(fields, first, count)
        !! Appends a data row, growing storage geometrically as needed.
        type(delimited_fields), intent(inout) :: fields
            !! The field collection.
        integer(int32), intent(in) :: first
            !! The index of the row's first field.
        integer(int32), intent(in) :: count
            !! The number of fields in the row.

        integer(int32), allocatable, dimension(:) :: temp
        integer(int32) :: n

        n = fields%row_count
        if (n == size(fields%row_first)) then
            allocate(temp(2 * n))
            temp(1:n) = fields%row_first
            call move_alloc(temp, fields%row_first)
            allocate(temp(2 * n))
            temp(1:n) = fields%row_size
            call move_alloc(temp, fields%row_size)
        end if
        n = n + 1
        fields%row_first(n) = first
        fields%row_size(n) = count
        fields%row_count = n
    end subroutine

! ------------------------------------------------------------------------------
    pure function cell_index(fields, row, col) result(rst)
        !! Locates the field holding the requested cell.
        type(delimited_fields), intent(in) :: fields
            !! The field collection.
        integer(int32), intent(in) :: row
            !! The data row index.
        integer(int32), intent(in) :: col
            !! The column index.
        integer(int32) :: rst
            !! The field index, or 0 if the row is too short to hold the cell.

        if (col > fields%row_size(row)) then
            rst = 0
        else
            rst = fields%row_first(row) + col - 1
        end if
    end function

! ------------------------------------------------------------------------------
    pure function get_field_text(buffer, fields, k) result(rst)
        !! Extracts a field's text, replacing escaped ("") quote pairs within
        !! quoted fields by single quotes.
        character(len = *), intent(in) :: buffer
            !! The file contents.
        type(delimited_fields), intent(in) :: fields
            !! The field collection.
        integer(int32), intent(in) :: k
            !! The field index.
        character(len = :), allocatable :: rst
            !! The field text.

        integer(int32) :: i, n

        rst = buffer(fields%first(k):fields%last(k))
        if (.not.fields%quoted(k) .or. index(rst, QUOTE // QUOTE) == 0) return
        n = 0
        i = 1
        do while (i <= len(rst))
            n = n + 1
            rst(n:n) = rst(i:i)
            if (rst(i:i) == QUOTE) i = i + 1
            i = i + 1
        end do
        rst = rst(1:n)
    end function

! ------------------------------------------------------------------------------
    subroutine add_inferred_column(table, buffer, fields, col, name)
        !! Infers a column's type from its contents and appends it to a table.
        type(data_table), intent(inout) :: table
            !! The table to which the column is appended.
        character(len = *), intent(in) :: buffer
            !! The file contents.
        type(delimited_fields), intent(in) :: fields
            !! The parsed data fields.
        integer(int32), intent(in) :: col
            !! The column index.
        character(len = *), intent(in) :: name
            !! The column name.

        real(real64), allocatable, dimension(:) :: values
        type(string), allocatable, dimension(:) :: texts
        type(integer_column) :: icol
        type(real_column) :: rcol
        type(logical_column) :: lcol
        type(string_column) :: scol
        integer(int32) :: i, k, nrows, n_logical, n_numeric
        logical :: all_integer, any_missing, is_text, empty
        real(real64) :: nan

        nrows = fields%row_count
        nan = ieee_value(nan, ieee_quiet_nan)
        allocate(values(nrows))
        n_logical = 0
        n_numeric = 0
        all_integer = .true.
        any_missing = .false.
        is_text = .false.

        ! Numeric and logical values are parsed only once; logical entries are
        ! held as 1 or 0 until the final type is known.
        do i = 1, nrows
            k = cell_index(fields, i, col)
            empty = k == 0
            if (.not.empty) empty = fields%last(k) < fields%first(k)
            if (empty) then
                any_missing = .true.
                values(i) = nan
                cycle
            end if
            select case (classify_field(buffer(fields%first(k):fields%last(k)), values(i)))
            case (FIELD_LOGICAL)
                n_logical = n_logical + 1
            case (FIELD_INTEGER)
                n_numeric = n_numeric + 1
            case (FIELD_REAL)
                n_numeric = n_numeric + 1
                all_integer = .false.
            case default
                is_text = .true.
                exit
            end select
        end do
        if (n_logical > 0 .and. (n_numeric > 0 .or. any_missing)) is_text = .true.

        if (is_text) then
            allocate(texts(nrows))
            do i = 1, nrows
                k = cell_index(fields, i, col)
                if (k == 0) then
                    texts(i) = ""
                else
                    texts(i) = get_field_text(buffer, fields, k)
                end if
            end do
            scol = string_column(texts)
            call scol%set_name(name)
            call table%add_column(scol)
        else if (n_logical > 0) then
            lcol = logical_column(values /= 0.0d0)
            call lcol%set_name(name)
            call table%add_column(lcol)
        else if (n_numeric > 0 .and. all_integer .and. .not.any_missing) then
            icol = integer_column(nint(values, int32))
            call icol%set_name(name)
            call table%add_column(icol)
        else
            rcol = real_column(values)
            call rcol%set_name(name)
            call table%add_column(rcol)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    function classify_field(text, value) result(rst)
        !! Determines whether a non-empty field is logical, integer, real, or
        !! text, and parses its value.
        character(len = *), intent(in) :: text
            !! The field text; must not be empty.
        real(real64), intent(out) :: value
            !! The parsed value: 1 or 0 for logical fields, the numeric value
            !! for numeric fields, and 0 for text.
        integer(int32) :: rst
            !! [[FIELD_LOGICAL]], [[FIELD_INTEGER]], [[FIELD_REAL]], or
            !! [[FIELD_TEXT]].

        character(len = 9) :: lower
        integer(int64) :: ival
        integer(int32) :: i, n, ndigits, nexp, ios
        logical :: is_integer

        n = len(text)
        value = 0.0d0
        rst = FIELD_TEXT

        if (n <= len(lower)) then
            lower = text
            do i = 1, n
                if (lge(lower(i:i), "A") .and. lle(lower(i:i), "Z")) &
                    lower(i:i) = achar(iachar(lower(i:i)) + 32)
            end do
            select case (lower)
            case ("t", "true", ".t.", ".true.")
                value = 1.0d0
                rst = FIELD_LOGICAL
                return
            case ("f", "false", ".f.", ".false.")
                rst = FIELD_LOGICAL
                return
            case ("nan", "+nan", "-nan")
                value = ieee_value(value, ieee_quiet_nan)
                rst = FIELD_REAL
                return
            case ("inf", "+inf", "infinity", "+infinity")
                value = ieee_value(value, ieee_positive_inf)
                rst = FIELD_REAL
                return
            case ("-inf", "-infinity")
                value = ieee_value(value, ieee_negative_inf)
                rst = FIELD_REAL
                return
            end select
        end if

        ! Validate [sign] digits [. digits] [(e|d) [sign] digits] before reading
        ! because list-directed input accepts many non-numeric forms.
        i = 1
        if (text(1:1) == "+" .or. text(1:1) == "-") i = 2
        ndigits = 0
        is_integer = .true.
        call skip_digits(text, i, ndigits)
        if (i <= n) then
            if (text(i:i) == ".") then
                is_integer = .false.
                i = i + 1
                call skip_digits(text, i, ndigits)
            end if
        end if
        if (ndigits == 0) return
        if (i <= n) then
            if (scan(text(i:i), "eEdD") == 0) return
            is_integer = .false.
            i = i + 1
            if (i <= n) then
                if (text(i:i) == "+" .or. text(i:i) == "-") i = i + 1
            end if
            nexp = 0
            call skip_digits(text, i, nexp)
            if (nexp == 0 .or. i <= n) return
        end if

        if (is_integer .and. ndigits <= 18) then
            read(text, *, iostat = ios) ival
            if (ios == 0 .and. ival >= -huge(1_int32) - 1 .and. &
                ival <= huge(1_int32)) &
            then
                value = real(ival, real64)
                rst = FIELD_INTEGER
                return
            end if
        end if
        read(text, *, iostat = ios) value
        if (ios == 0) then
            rst = FIELD_REAL
        else
            value = 0.0d0
        end if
    end function

! ------------------------------------------------------------------------------
    pure subroutine skip_digits(text, i, count)
        !! Advances past a run of decimal digits.
        character(len = *), intent(in) :: text
            !! The text being scanned.
        integer(int32), intent(inout) :: i
            !! The scan position.
        integer(int32), intent(inout) :: count
            !! Incremented by the number of digits skipped.

        do while (i <= len(text))
            if (llt(text(i:i), "0") .or. lgt(text(i:i), "9")) exit
            i = i + 1
            count = count + 1
        end do
    end subroutine

! ******************************************************************************
! WRITING
! ------------------------------------------------------------------------------
    subroutine write_delimited_table(fname, table, delimiter, write_header, &
        real_fmt)
        !! Writes a [[data_table]] to a delimited text file, replacing any
        !! existing file.
        !!
        !! The file contents are assembled in memory and written to disk with
        !! a single operation.  Rows end with a line feed.  NaN values in real
        !! columns are written as empty (missing) fields.  Logical values are
        !! written as T or F.  Strings and column names are enclosed in double
        !! quotes, with embedded quotes doubled, when they contain the
        !! delimiter, a quote, or a line break, or begin or end with a space
        !! or tab.  The output can be read with [[read_delimited_table]];
        !! note, however, that string columns whose values all resemble
        !! numbers or logicals will be read back with the inferred type.
        !!
        !! The following errors terminate the program via error stop:
        !!
        !! - [[FS_INVALID_INPUT_ERROR]]: the delimiter is invalid, or the file
        !!   cannot be opened or written.
        character(len = *), intent(in) :: fname
            !! The name of the file to write.
        class(data_table), intent(in) :: table
            !! The [[data_table]] to write.
        character, intent(in), optional :: delimiter
            !! The field delimiter.  The default is a comma.  The delimiter
            !! must not be a double quote, carriage return, or line feed.
        logical, intent(in), optional :: write_header
            !! Set to true to write the column names as the first row, or
            !! false to omit them.  The default is true.
        character(len = *), intent(in), optional :: real_fmt
            !! An optional Fortran edit descriptor, without parentheses, used
            !! for real values (e.g. "ES24.16E3" for full double precision).
            !! If omitted, the [[real_column]] default formatting is used.

        character :: delim
        logical :: header
        character(len = :), allocatable :: buffer
        type(column_pointer), allocatable, dimension(:) :: cols
        integer(int32) :: i, j, ncols, nrows, length, unit, ios

        delim = ","
        if (present(delimiter)) delim = delimiter
        if (delim == QUOTE .or. delim == CR .or. delim == LF) &
            error stop FS_INVALID_INPUT_ERROR
        header = .true.
        if (present(write_header)) header = write_header

        ncols = table%column_count()
        nrows = table%row_count()
        allocate(cols(ncols))
        do j = 1, ncols
            cols(j)%ptr => table%get_column(j)
        end do

        allocate(character(len = 4096) :: buffer)
        length = 0
        if (header .and. ncols > 0) then
            do j = 1, ncols
                if (j > 1) call append_text(buffer, length, delim)
                if (allocated(cols(j)%ptr%m_name)) &
                    call append_text(buffer, length, &
                        quote_field(cols(j)%ptr%m_name, delim))
            end do
            call append_text(buffer, length, LF)
        end if

        do i = 1, nrows
            do j = 1, ncols
                if (j > 1) call append_text(buffer, length, delim)
                select type (col => cols(j)%ptr)
                type is (real_column)
                    if (ieee_is_nan(col%get(i))) cycle
                    if (present(real_fmt)) then
                        call append_text(buffer, length, col%to_string(i, real_fmt))
                    else
                        call append_text(buffer, length, col%to_string(i))
                    end if
                type is (string_column)
                    call append_text(buffer, length, quote_field(col%get(i), delim))
                class default
                    call append_text(buffer, length, col%to_string(i))
                end select
            end do
            call append_text(buffer, length, LF)
        end do

        open(newunit = unit, file = fname, access = "stream", &
            form = "unformatted", action = "write", status = "replace", &
            iostat = ios)
        if (ios /= 0) error stop FS_INVALID_INPUT_ERROR
        if (length > 0) write(unit, iostat = ios) buffer(1:length)
        close(unit)
        if (ios /= 0) error stop FS_INVALID_INPUT_ERROR
    end subroutine

! ------------------------------------------------------------------------------
    pure subroutine append_text(buffer, length, text)
        !! Appends text to a buffer, growing it geometrically as needed.
        character(len = :), allocatable, intent(inout) :: buffer
            !! The buffer.
        integer(int32), intent(inout) :: length
            !! The number of characters in use.
        character(len = *), intent(in) :: text
            !! The text to append.

        character(len = :), allocatable :: temp
        integer(int32) :: n

        n = len(text)
        if (length + n > len(buffer)) then
            allocate(character(len = max(2 * len(buffer), length + n)) :: temp)
            temp(1:length) = buffer(1:length)
            call move_alloc(temp, buffer)
        end if
        buffer(length + 1:length + n) = text
        length = length + n
    end subroutine

! ------------------------------------------------------------------------------
    pure function quote_field(text, delim) result(rst)
        !! Encloses text in double quotes, doubling embedded quotes, if it
        !! could not otherwise be read back verbatim.
        character(len = *), intent(in) :: text
            !! The text to write.
        character, intent(in) :: delim
            !! The field delimiter.
        character(len = :), allocatable :: rst
            !! The text, quoted if necessary.

        integer(int32) :: i, n
        logical :: needs_quotes

        n = len(text)
        needs_quotes = .false.
        if (n > 0) needs_quotes = scan(text, delim // QUOTE // CR // LF) > 0 &
            .or. text(1:1) == " " .or. text(1:1) == TAB &
            .or. text(n:n) == " " .or. text(n:n) == TAB
        if (.not.needs_quotes) then
            rst = text
            return
        end if

        rst = QUOTE
        do i = 1, n
            if (text(i:i) == QUOTE) then
                rst = rst // QUOTE // QUOTE
            else
                rst = rst // text(i:i)
            end if
        end do
        rst = rst // QUOTE
    end function

! ------------------------------------------------------------------------------
end module
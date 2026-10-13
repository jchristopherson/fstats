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

module fstats_data_columns
    use iso_fortran_env
    use fstats_errors
    use strings
    use, intrinsic :: ieee_arithmetic
    implicit none
    private
    public :: FS_COLUMN_INTEGER
    public :: FS_COLUMN_REAL
    public :: FS_COLUMN_LOGICAL
    public :: FS_COLUMN_STRING
    public :: data_column
    public :: integer_column
    public :: real_column
    public :: logical_column
    public :: string_column
    public :: size
    public :: string
    public :: ieee_is_finite, ieee_is_nan, ieee_is_negative, ieee_is_normal
    public :: abs, maxval, minval, sin, cos, tan, asin, acos, atan
    public :: exp, log, log10, sum, product, sqrt
    public :: operator(+), operator(-), operator(*), operator(/)

    integer(int32), parameter :: FS_COLUMN_INTEGER = 1
        !! Denotes a column of integer values.
    integer(int32), parameter :: FS_COLUMN_REAL = 2
        !! Denotes a column of real (floating-point) values.
    integer(int32), parameter :: FS_COLUMN_LOGICAL = 3
        !! Denotes a column of logical values.
    integer(int32), parameter :: FS_COLUMN_STRING = 4
        !! Denotes a column of text strings.

    type, abstract :: data_column
        !! Defines a column of data within a [[data_table]] object.
        character(len = :), allocatable :: m_name
            !! A header, or column name.
    contains
        procedure(data_column_integer_query), deferred, public, pass :: count
            !! Gets the number of items stored in the column.
        procedure(data_column_integer_query), deferred, public, pass :: get_data_type
            !! Gets the data type stored by the column:
            !!
            !! - [[FS_COLUMN_INTEGER]]
            !!
            !! - [[FS_COLUMN_REAL]]
            !!
            !! - [[FS_COLUMN_LOGICAL]]
            !!
            !! - [[FS_COLUMN_STRING]]
        procedure(data_column_array_copy), deferred, public, pass :: copy
            !! Copies the contents of the column to an array.
        procedure(data_column_array_view), deferred, public, pass :: view
            !! Borrows a mutable array pointer without copying the column.
            !! The actual column must have TARGET or be an associated pointer.
            !! The caller must not deallocate the returned pointer.
            !! Reallocation of the array or destruction of the column invalidates
            !! the view; discard it before either operation and reacquire afterward.
        procedure(data_column_slice), deferred, public, pass :: slice
            !! Borrows a mutable array pointer to a slice of the column without
            !! copying the column.  The actual column must have TARGET or be an 
            !! associated pointer.  The caller must not deallocate the returned 
            !! pointer.  Reallocation of the array or destruction of the column 
            !! invalidates the view; discard it before either operation and 
            !! reacquire afterward.
        procedure(data_column_to_string), deferred, public, pass :: to_string
            !! Converts the requested value to a string.
        procedure(data_column_append), deferred, public, pass :: append
            !! Appends one value, initializing an unallocated column if needed.
            !! Discard existing views before this operation.
        procedure(data_column_insert), deferred, public, pass :: insert
            !! Inserts one value before row, in the range 1 through count() + 1.
            !! Discard existing views before this operation.
        procedure(data_column_remove), deferred, public, pass :: remove
            !! Removes one row, in the range 1 through count().
            !! Discard existing views before this operation.
        procedure, public :: get_name => dc_get_name
        procedure, public :: set_name => dc_set_name
    end type

    abstract interface
        pure function data_column_integer_query(this) result(rst)
            !! Defines an integer property query on a [[data_column]] object.
            use iso_fortran_env, only : int32
            import data_column
            class(data_column), intent(in) :: this
                !! The [[data_column]] object.
            integer(int32) :: rst
                !! The property value.
        end function

        pure subroutine data_column_array_copy(this, dest)
            !! Defines a mechanism to retrieve a copy of an array.
            !! The destination must match the column's type and count.
            !! An unallocated column is treated as empty.
            import data_column
            class(data_column), intent(in) :: this
                !! The [[data_column]] object.
            class(*), intent(inout), dimension(:) :: dest
                !! The destination array, overwritten in full on success.
        end subroutine

        function data_column_array_view(this) result(rst)
            !! Returns a borrowed pointer to the column's array.
            import data_column
            class(data_column), target, intent(inout) :: this
                !! The [[data_column]] object that owns the array.
            class(*), pointer, dimension(:) :: rst
                !! The mutable view, or null if the array is unallocated.
        end function

        function data_column_to_string(this, row, fmt) result(rst)
            !! Returns the value in the requested row as a string.
            use iso_fortran_env, only : int32
            import data_column
            class(data_column), intent(in) :: this
                !! The [[data_column]] object.
            integer(int32), intent(in) :: row
                !! The row index.
            character(len = *), intent(in), optional :: fmt
                !! The formatting string.
            character(len = :), allocatable :: rst
                !! The resulting string.
        end function
        subroutine data_column_append(this, x)
            !! Appends a scalar of the column's type; strings also accept character.
            import data_column
            class(data_column), intent(inout) :: this
            class(*), intent(in) :: x
        end subroutine

        subroutine data_column_insert(this, row, x)
            !! Inserts a scalar of the column's type; strings also accept character.
            import data_column, int32
            class(data_column), intent(inout) :: this
            integer(int32), intent(in) :: row
            class(*), intent(in) :: x
        end subroutine

        subroutine data_column_remove(this, row)
            import data_column, int32
            class(data_column), intent(inout) :: this
            integer(int32), intent(in) :: row
        end subroutine

        ! ------------------------------------------------------------------------------
        function data_column_slice(this, first, last) result(rst)
            !! Returns a borrowed, mutable slice of this column's array.
            use iso_fortran_env, only : int32
            import data_column
            class(data_column), target, intent(inout) :: this
                !! The [[data_column]] object.
            integer(int32), intent(in) :: first
                !! The one-based starting index of the slice.
            integer(int32), intent(in) :: last
                !! The one-based ending index of the slice.
            class(*), pointer, dimension(:) :: rst
                !! The array pointer, or null if unallocated.
        end function
    end interface

    type, extends(data_column) :: integer_column
        !! Defines a column of integer values.
        integer(int32), private, allocatable, dimension(:) :: m_values
    contains
        procedure, public :: count => ic_count
        procedure, public :: get_data_type => ic_get_type
        procedure, public :: get => ic_get
        procedure, public :: set => ic_set
        procedure, public :: copy => ic_copy
        procedure, public :: view => ic_view
        procedure, public :: slice => ic_slice
        procedure, public :: to_string => ic_to_string
        procedure, public :: append => ic_append
        procedure, public :: insert => ic_insert
        procedure, public :: remove => ic_remove
    end type

    interface integer_column
        module procedure :: ic_init_empty
        module procedure :: ic_init_array
    end interface

    type, extends(data_column) :: real_column
        !! Defines a column of double-precision, floating-point values.
        real(real64), private, allocatable, dimension(:) :: m_values
    contains
        procedure, public :: count => rc_count
        procedure, public :: get_data_type => rc_get_type
        procedure, public :: get => rc_get
        procedure, public :: set => rc_set
        procedure, public :: copy => rc_copy
        procedure, public :: view => rc_view
        procedure, public :: slice => rc_slice
        procedure, public :: to_string => rc_to_string
        procedure, public :: append => rc_append
        procedure, public :: insert => rc_insert
        procedure, public :: remove => rc_remove
    end type

    interface real_column
        module procedure :: rc_init_empty
        module procedure :: rc_init_array
    end interface

    type, extends(data_column) :: logical_column
        !! Defines a column of logical values.
        logical, private, allocatable, dimension(:) :: m_values
    contains
        procedure, public :: count => lc_count
        procedure, public :: get_data_type => lc_get_type
        procedure, public :: get => lc_get
        procedure, public :: set => lc_set
        procedure, public :: copy => lc_copy
        procedure, public :: view => lc_view
        procedure, public :: slice => lc_slice
        procedure, public :: to_string => lc_to_string
        procedure, public :: append => lc_append
        procedure, public :: insert => lc_insert
        procedure, public :: remove => lc_remove
    end type

    interface logical_column
        module procedure :: lc_init_empty
        module procedure :: lc_init_array
    end interface

    type, extends(data_column) :: string_column
        !! Defines a column of string values.
        type(string), private, allocatable, dimension(:) :: m_values
    contains
        procedure, public :: count => sc_count
        procedure, public :: get_data_type => sc_get_type
        procedure, public :: get => sc_get
        procedure, public :: set => sc_set
        procedure, public :: copy => sc_copy
        procedure, public :: view => sc_view
        procedure, public :: slice => sc_slice
        procedure, public :: to_string => sc_to_string
        procedure, public :: append => sc_append
        procedure, public :: insert => sc_insert
        procedure, public :: remove => sc_remove
    end type

    interface string_column
        module procedure :: sc_init_empty
        module procedure :: sc_init_array
    end interface

! ------------------------------------------------------------------------------
    interface size
        module procedure :: dc_size
    end interface

    interface ieee_is_finite
        module procedure :: rc_ieee_is_finite
    end interface
    interface abs
        module procedure :: rc_abs
    end interface
    interface maxval
        module procedure :: rc_maxval
    end interface
    interface minval
        module procedure :: rc_minval
    end interface
    interface sin
        module procedure :: rc_sin
    end interface
    interface cos
        module procedure :: rc_cos
    end interface
    interface tan
        module procedure :: rc_tan
    end interface
    interface asin
        module procedure :: rc_asin
    end interface
    interface acos
        module procedure :: rc_acos
    end interface
    interface atan
        module procedure :: rc_atan
    end interface
    interface exp
        module procedure :: rc_exp
    end interface
    interface log
        module procedure :: rc_log
    end interface
    interface log10
        module procedure :: rc_log10
    end interface
    interface sum
        module procedure :: rc_sum
    end interface
    interface product
        module procedure :: rc_product
    end interface
    interface sqrt
        module procedure :: rc_sqrt
    end interface
    interface operator(+)
        module procedure :: rc_add, rc_add_scalar, scalar_add_rc, rc_positive
    end interface
    interface operator(-)
        module procedure :: rc_subtract, rc_subtract_scalar, scalar_subtract_rc, rc_negative
    end interface
    interface operator(*)
        module procedure :: rc_multiply, rc_multiply_scalar, scalar_multiply_rc
    end interface
    interface operator(/)
        module procedure :: rc_divide, rc_divide_scalar, scalar_divide_rc
    end interface
    interface ieee_is_nan
        module procedure :: rc_ieee_is_nan
    end interface
    interface ieee_is_negative
        module procedure :: rc_ieee_is_negative
    end interface
    interface ieee_is_normal
        module procedure :: rc_ieee_is_normal
    end interface

contains
! ******************************************************************************
! DATA_COLUMN MEMBERS
! ------------------------------------------------------------------------------
    pure function dc_get_name(this) result(rst)
        !! Gets the name associated with the column.
        class(data_column), intent(in) :: this
            !! The [[data_column]] object.
        character(len = :), allocatable :: rst
            !! The name.
        rst = this%m_name
    end function

! ------------------------------------------------------------------------------
    subroutine dc_set_name(this, x)
        !! Sets the name associated with the column.
        class(data_column), intent(inout) :: this
            !! The [[data_column]] object.
        character(len = *), intent(in) :: x
            !! The name.
        this%m_name = x
    end subroutine

! ******************************************************************************
! INTEGER_COLUMN MEMBERS
! ------------------------------------------------------------------------------
    pure function ic_init_empty(n) result(rst)
        !! Constructs an empty [[integer_column]] of the requested size.
        integer(int32), intent(in) :: n
            !! The number of rows to create.
        type(integer_column) :: rst
            !! The newly constructed [[integer_column]] object.
        
        if (n < 0) error stop FS_INVALID_INPUT_ERROR
        allocate(rst%m_values(n), source = 0)
    end function

! ------------------------------------------------------------------------------
    pure function ic_init_array(x) result(rst)
        !! Constructs a new [[integer_column]] by copying the contents of the
        !! supplied array.
        integer(int32), intent(in), dimension(:) :: x
            !! The array to copy.
        type(integer_column) :: rst
            !! The newly constructed [[integer_column]] object.

        allocate(rst%m_values(size(x)), source = x)
    end function

! ------------------------------------------------------------------------------
    pure function ic_count(this) result(rst)
        !! Gets the number of values stored in this column.
        class(integer_column), intent(in) :: this
            !! The [[integer_column]] object.
        integer(int32) :: rst
            !! The number of values.

        if (allocated(this%m_values)) then
            rst = size(this%m_values)
        else
            rst = 0
        end if
    end function

! ------------------------------------------------------------------------------
    pure function ic_get_type(this) result(rst)
        !! Gets the type stored by this column.
        class(integer_column), intent(in) :: this
            !! The [[integer_column]] object.
        integer(int32) :: rst
            !! The type stored by this column.
        rst = FS_COLUMN_INTEGER
    end function

! ------------------------------------------------------------------------------
    pure function ic_get(this, row) result(rst)
        !! Gets the requested value from the column.
        class(integer_column), intent(in) :: this
            !! The [[integer_column]] object.
        integer(int32), intent(in) :: row
            !! The row index of the value to retrieve.
        integer(int32) :: rst
            !! The value stored at the requested row.

        if (.not.allocated(this%m_values)) error stop FS_UNINITIALIZED_OBJECT_ERROR
        if (row < 1 .or. row > this%count()) error stop FS_INVALID_INPUT_ERROR
        rst = this%m_values(row)
    end function

! --------------------
    subroutine ic_set(this, row, x)
        !! Sets the requested value into the column.
        class(integer_column), intent(inout) :: this
            !! The [[integer_column]] object.
        integer(int32), intent(in) :: row
            !! The row index.
        integer(int32), intent(in) :: x
            !! The value stored at the requested row.

        if (.not.allocated(this%m_values)) error stop FS_UNINITIALIZED_OBJECT_ERROR
        if (row < 1 .or. row > this%count()) error stop FS_INVALID_INPUT_ERROR
        this%m_values(row) = x
    end subroutine

! ------------------------------------------------------------------------------
    pure subroutine ic_copy(this, dest)
        !! Copies the contents of this column to the supplied destination array.
        class(integer_column), intent(in) :: this
            !! The [[integer_column]] object.
        class(*), intent(inout), dimension(:) :: dest
            !! The destination array.

        if (size(dest) /= this%count()) error stop FS_ARRAY_SIZE_ERROR
        select type (dest)
        type is (integer(int32))
            if (allocated(this%m_values)) dest = this%m_values
        class default
            error stop FS_INVALID_ARGUMENT_ERROR
        end select
    end subroutine

! ------------------------------------------------------------------------------
    function ic_view(this) result(rst)
        !! Returns a borrowed, mutable view of this column's array.
        class(integer_column), target, intent(inout) :: this
            !! The [[integer_column]] object.
        class(*), pointer, dimension(:) :: rst
            !! The array pointer, or null if unallocated.

        nullify(rst)
        if (allocated(this%m_values)) rst => this%m_values
    end function

! ------------------------------------------------------------------------------
    function ic_slice(this, first, last) result(rst)
        !! Returns a borrowed, mutable slice of this column's array.
        class(integer_column), target, intent(inout) :: this
            !! The [[integer_column]] object.
        integer(int32), intent(in) :: first
            !! The one-based starting index of the slice.
        integer(int32), intent(in) :: last
            !! The one-based ending index of the slice.
        class(*), pointer, dimension(:) :: rst
            !! The array pointer, or null if unallocated.

        nullify(rst)
        if (first > last) error stop FS_INVALID_ARGUMENT_ERROR
        if (first < 1) error stop FS_INVALID_INPUT_ERROR
        if (last > this%count()) error stop FS_INVALID_INPUT_ERROR
        if (allocated(this%m_values)) rst => this%m_values(first:last)
    end function

! ------------------------------------------------------------------------------
    function ic_to_string(this, row, fmt) result(rst)
        !! Returns the value in the requested row as a string.
        class(integer_column), intent(in) :: this
            !! The [[integer_column]] object.
        integer(int32), intent(in) :: row
            !! The row index.
        character(len = *), intent(in), optional :: fmt
            !! The formatting string (e.g. "I0").
        character(len = :), allocatable :: rst
            !! The resulting string.

        character(len = :), allocatable :: frmt, buffer
        
        if (.not.allocated(this%m_values)) error stop FS_UNINITIALIZED_OBJECT_ERROR
        if (row < 1 .or. row > this%count()) error stop FS_INVALID_INPUT_ERROR

        allocate(character(len = 50) :: buffer)
        if (present(fmt)) then
            frmt = "(" // fmt // ")"
        else
            frmt = "(I0)"
        end if
        write(buffer, frmt) this%m_values(row)
        rst = trim(adjustl(buffer))
    end function

! ------------------------------------------------------------------------------
    subroutine ic_append(this, x)
        class(integer_column), intent(inout) :: this
        class(*), intent(in) :: x

        call this%insert(this%count() + 1, x)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine ic_insert(this, row, x)
        class(integer_column), intent(inout) :: this
        integer(int32), intent(in) :: row
        class(*), intent(in) :: x
        integer(int32), allocatable, dimension(:) :: values
        integer(int32) :: n

        n = this%count()
        if (row < 1 .or. row > n + 1) error stop FS_INVALID_INPUT_ERROR
        select type (x)
        type is (integer(int32))
            allocate(values(n + 1))
            if (allocated(this%m_values)) then
                values(:row - 1) = this%m_values(:row - 1)
                values(row + 1:) = this%m_values(row:)
            end if
            values(row) = x
        class default
            error stop FS_INVALID_ARGUMENT_ERROR
        end select
        call move_alloc(values, this%m_values)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine ic_remove(this, row)
        class(integer_column), intent(inout) :: this
        integer(int32), intent(in) :: row
        integer(int32), allocatable, dimension(:) :: values
        integer(int32) :: n

        n = this%count()
        if (row < 1 .or. row > n) error stop FS_INVALID_INPUT_ERROR
        allocate(values(n - 1))
        values(:row - 1) = this%m_values(:row - 1)
        values(row:) = this%m_values(row + 1:)
        call move_alloc(values, this%m_values)
    end subroutine

! ******************************************************************************
! REAL_COLUMN MEMBERS
! ------------------------------------------------------------------------------
    pure function rc_init_empty(n) result(rst)
        !! Constructs an empty [[real_column]] of the requested size.
        integer(int32), intent(in) :: n
            !! The number of rows to create.
        type(real_column) :: rst
            !! The newly constructed [[real_column]] object.
        
        if (n < 0) error stop FS_INVALID_INPUT_ERROR
        allocate(rst%m_values(n), source = 0.0d0)
    end function

! ------------------------------------------------------------------------------
    pure function rc_init_array(x) result(rst)
        !! Constructs a new [[real_column]] by copying the contents of the
        !! supplied array.
        real(real64), intent(in), dimension(:) :: x
            !! The array to copy.
        type(real_column) :: rst
            !! The newly constructed [[real_column]] object.

        allocate(rst%m_values(size(x)), source = x)
    end function

! ------------------------------------------------------------------------------
    pure function rc_count(this) result(rst)
        !! Gets the number of values stored in this column.
        class(real_column), intent(in) :: this
            !! The [[real_column]] object.
        integer(int32) :: rst
            !! The number of values.

        if (allocated(this%m_values)) then
            rst = size(this%m_values)
        else
            rst = 0
        end if
    end function

! ------------------------------------------------------------------------------
    pure function rc_get_type(this) result(rst)
        !! Gets the type stored by this column.
        class(real_column), intent(in) :: this
            !! The [[real_column]] object.
        integer(int32) :: rst
            !! The type stored by this column.
        rst = FS_COLUMN_REAL
    end function

! ------------------------------------------------------------------------------
    pure function rc_get(this, row) result(rst)
        !! Gets the requested value from the column.
        class(real_column), intent(in) :: this
            !! The [[real_column]] object.
        integer(int32), intent(in) :: row
            !! The row index of the value to retrieve.
        real(real64) :: rst
            !! The value stored at the requested row.

        if (.not.allocated(this%m_values)) error stop FS_UNINITIALIZED_OBJECT_ERROR
        if (row < 1 .or. row > this%count()) error stop FS_INVALID_INPUT_ERROR
        rst = this%m_values(row)
    end function

! --------------------
    subroutine rc_set(this, row, x)
        !! Sets the requested value into the column.
        class(real_column), intent(inout) :: this
            !! The [[real_column]] object.
        integer(int32), intent(in) :: row
            !! The row index.
        real(real64), intent(in) :: x
            !! The value stored at the requested row.

        if (.not.allocated(this%m_values)) error stop FS_UNINITIALIZED_OBJECT_ERROR
        if (row < 1 .or. row > this%count()) error stop FS_INVALID_INPUT_ERROR
        this%m_values(row) = x
    end subroutine

! ------------------------------------------------------------------------------
    pure subroutine rc_copy(this, dest)
        !! Copies the contents of this column to the supplied destination array.
        class(real_column), intent(in) :: this
            !! The [[real_column]] object.
        class(*), intent(inout), dimension(:) :: dest
            !! The destination array.

        if (size(dest) /= this%count()) error stop FS_ARRAY_SIZE_ERROR
        select type (dest)
        type is (real(real64))
            if (allocated(this%m_values)) dest = this%m_values
        class default
            error stop FS_INVALID_ARGUMENT_ERROR
        end select
    end subroutine

! ------------------------------------------------------------------------------
    function rc_view(this) result(rst)
        !! Returns a borrowed, mutable view of this column's array.
        class(real_column), target, intent(inout) :: this
            !! The [[real_column]] object.
        class(*), pointer, dimension(:) :: rst
            !! The array pointer, or null if unallocated.

        nullify(rst)
        if (allocated(this%m_values)) rst => this%m_values
    end function

! ------------------------------------------------------------------------------
    function rc_slice(this, first, last) result(rst)
        !! Returns a borrowed, mutable slice of this column's array.
        class(real_column), target, intent(inout) :: this
            !! The [[real_column]] object.
        integer(int32), intent(in) :: first
            !! The one-based starting index of the slice.
        integer(int32), intent(in) :: last
            !! The one-based ending index of the slice.
        class(*), pointer, dimension(:) :: rst
            !! The array pointer, or null if unallocated.

        nullify(rst)
        if (first > last) error stop FS_INVALID_ARGUMENT_ERROR
        if (first < 1) error stop FS_INVALID_INPUT_ERROR
        if (last > this%count()) error stop FS_INVALID_INPUT_ERROR
        if (allocated(this%m_values)) rst => this%m_values(first:last)
    end function

! ------------------------------------------------------------------------------
    function rc_to_string(this, row, fmt) result(rst)
        !! Returns the value in the requested row as a string.
        class(real_column), intent(in) :: this
            !! The [[real_column]] object.
        integer(int32), intent(in) :: row
            !! The row index.
        character(len = *), intent(in), optional :: fmt
            !! The formatting string (e.g. "F0.6").
        character(len = :), allocatable :: rst
            !! The resulting string.

        character(len = :), allocatable :: frmt, buffer
        real(real64) :: val
        
        if (.not.allocated(this%m_values)) error stop FS_UNINITIALIZED_OBJECT_ERROR
        if (row < 1 .or. row > this%count()) error stop FS_INVALID_INPUT_ERROR

        allocate(character(len = 50) :: buffer)
        val = this%m_values(row)
        if (present(fmt)) then
            frmt = "(" // fmt // ")"
        else
            if (abs(val) < 1.0d-3 .or. abs(val) > 1.0d3) then
                frmt = "(ES18.8E3)"
            else
                frmt = "(F18.8)"
            end if
        end if
        write(buffer, frmt) val
        rst = trim_trailing_zeros(trim(adjustl(buffer)))
    end function

! ------------------------------------------------------------------------------
    subroutine rc_append(this, x)
        class(real_column), intent(inout) :: this
        class(*), intent(in) :: x

        call this%insert(this%count() + 1, x)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine rc_insert(this, row, x)
        class(real_column), intent(inout) :: this
        integer(int32), intent(in) :: row
        class(*), intent(in) :: x
        real(real64), allocatable, dimension(:) :: values
        integer(int32) :: n

        n = this%count()
        if (row < 1 .or. row > n + 1) error stop FS_INVALID_INPUT_ERROR
        select type (x)
        type is (real(real64))
            allocate(values(n + 1))
            if (allocated(this%m_values)) then
                values(:row - 1) = this%m_values(:row - 1)
                values(row + 1:) = this%m_values(row:)
            end if
            values(row) = x
        class default
            error stop FS_INVALID_ARGUMENT_ERROR
        end select
        call move_alloc(values, this%m_values)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine rc_remove(this, row)
        class(real_column), intent(inout) :: this
        integer(int32), intent(in) :: row
        real(real64), allocatable, dimension(:) :: values
        integer(int32) :: n

        n = this%count()
        if (row < 1 .or. row > n) error stop FS_INVALID_INPUT_ERROR
        allocate(values(n - 1))
        values(:row - 1) = this%m_values(:row - 1)
        values(row:) = this%m_values(row + 1:)
        call move_alloc(values, this%m_values)
    end subroutine

! ******************************************************************************
! LOGICAL_COLUMN MEMBERS
! ------------------------------------------------------------------------------
    pure function lc_init_empty(n) result(rst)
        !! Constructs an empty [[logical_column]] of the requested size.
        integer(int32), intent(in) :: n
            !! The number of rows to create.
        type(logical_column) :: rst
            !! The newly constructed [[logical_column]] object.
        
        if (n < 0) error stop FS_INVALID_INPUT_ERROR
        allocate(rst%m_values(n), source = .false.)
    end function

! ------------------------------------------------------------------------------
    pure function lc_init_array(x) result(rst)
        !! Constructs a new [[logical_column]] by copying the contents of the
        !! supplied array.
        logical, intent(in), dimension(:) :: x
            !! The array to copy.
        type(logical_column) :: rst
            !! The newly constructed [[logical_column]] object.

        allocate(rst%m_values(size(x)), source = x)
    end function

! ------------------------------------------------------------------------------
    pure function lc_count(this) result(rst)
        !! Gets the number of values stored in this column.
        class(logical_column), intent(in) :: this
            !! The [[logical_column]] object.
        integer(int32) :: rst
            !! The number of values.

        if (allocated(this%m_values)) then
            rst = size(this%m_values)
        else
            rst = 0
        end if
    end function

! ------------------------------------------------------------------------------
    pure function lc_get_type(this) result(rst)
        !! Gets the type stored by this column.
        class(logical_column), intent(in) :: this
            !! The [[logical_column]] object.
        integer(int32) :: rst
            !! The type stored by this column.
        rst = FS_COLUMN_LOGICAL
    end function

! ------------------------------------------------------------------------------
    pure function lc_get(this, row) result(rst)
        !! Gets the requested value from the column.
        class(logical_column), intent(in) :: this
            !! The [[logical_column]] object.
        integer(int32), intent(in) :: row
            !! The row index of the value to retrieve.
        logical :: rst
            !! The value stored at the requested row.

        if (.not.allocated(this%m_values)) error stop FS_UNINITIALIZED_OBJECT_ERROR
        if (row < 1 .or. row > this%count()) error stop FS_INVALID_INPUT_ERROR
        rst = this%m_values(row)
    end function

! --------------------
    subroutine lc_set(this, row, x)
        !! Sets the requested value into the column.
        class(logical_column), intent(inout) :: this
            !! The [[logical_column]] object.
        integer(int32), intent(in) :: row
            !! The row index.
        logical, intent(in) :: x
            !! The value stored at the requested row.

        if (.not.allocated(this%m_values)) error stop FS_UNINITIALIZED_OBJECT_ERROR
        if (row < 1 .or. row > this%count()) error stop FS_INVALID_INPUT_ERROR
        this%m_values(row) = x
    end subroutine

! ------------------------------------------------------------------------------
    pure subroutine lc_copy(this, dest)
        !! Copies the contents of this column to the supplied destination array.
        class(logical_column), intent(in) :: this
            !! The [[logical_column]] object.
        class(*), intent(inout), dimension(:) :: dest
            !! The destination array.

        if (size(dest) /= this%count()) error stop FS_ARRAY_SIZE_ERROR
        select type (dest)
        type is (logical)
            if (allocated(this%m_values)) dest = this%m_values
        class default
            error stop FS_INVALID_ARGUMENT_ERROR
        end select
    end subroutine

! ------------------------------------------------------------------------------
    function lc_view(this) result(rst)
        !! Returns a borrowed, mutable view of this column's array.
        class(logical_column), target, intent(inout) :: this
            !! The [[logical_column]] object.
        class(*), pointer, dimension(:) :: rst
            !! The array pointer, or null if unallocated.

        nullify(rst)
        if (allocated(this%m_values)) rst => this%m_values
    end function

! ------------------------------------------------------------------------------
    function lc_slice(this, first, last) result(rst)
        !! Returns a borrowed, mutable slice of this column's array.
        class(logical_column), target, intent(inout) :: this
            !! The [[logical_column]] object.
        integer(int32), intent(in) :: first
            !! The one-based starting index of the slice.
        integer(int32), intent(in) :: last
            !! The one-based ending index of the slice.
        class(*), pointer, dimension(:) :: rst
            !! The array pointer, or null if unallocated.

        nullify(rst)
        if (first > last) error stop FS_INVALID_ARGUMENT_ERROR
        if (first < 1) error stop FS_INVALID_INPUT_ERROR
        if (last > this%count()) error stop FS_INVALID_INPUT_ERROR
        if (allocated(this%m_values)) rst => this%m_values(first:last)
    end function

! ------------------------------------------------------------------------------
    function lc_to_string(this, row, fmt) result(rst)
        !! Returns the value in the requested row as a string.
        class(logical_column), intent(in) :: this
            !! The [[logical_column]] object.
        integer(int32), intent(in) :: row
            !! The row index.
        character(len = *), intent(in), optional :: fmt
            !! The formatting string (e.g. "L1").
        character(len = :), allocatable :: rst
            !! The resulting string.

        character(len = :), allocatable :: frmt, buffer
        
        if (.not.allocated(this%m_values)) error stop FS_UNINITIALIZED_OBJECT_ERROR
        if (row < 1 .or. row > this%count()) error stop FS_INVALID_INPUT_ERROR

        allocate(character(len = 10) :: buffer)
        if (present(fmt)) then
            frmt = "(" // fmt // ")"
        else
            frmt = "(L1)"
        end if
        write(buffer, frmt) this%m_values(row)
        rst = trim(adjustl(buffer))
    end function

! ------------------------------------------------------------------------------
    subroutine lc_append(this, x)
        class(logical_column), intent(inout) :: this
        class(*), intent(in) :: x

        call this%insert(this%count() + 1, x)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine lc_insert(this, row, x)
        class(logical_column), intent(inout) :: this
        integer(int32), intent(in) :: row
        class(*), intent(in) :: x
        logical, allocatable, dimension(:) :: values
        integer(int32) :: n

        n = this%count()
        if (row < 1 .or. row > n + 1) error stop FS_INVALID_INPUT_ERROR
        select type (x)
        type is (logical)
            allocate(values(n + 1))
            if (allocated(this%m_values)) then
                values(:row - 1) = this%m_values(:row - 1)
                values(row + 1:) = this%m_values(row:)
            end if
            values(row) = x
        class default
            error stop FS_INVALID_ARGUMENT_ERROR
        end select
        call move_alloc(values, this%m_values)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine lc_remove(this, row)
        class(logical_column), intent(inout) :: this
        integer(int32), intent(in) :: row
        logical, allocatable, dimension(:) :: values
        integer(int32) :: n

        n = this%count()
        if (row < 1 .or. row > n) error stop FS_INVALID_INPUT_ERROR
        allocate(values(n - 1))
        values(:row - 1) = this%m_values(:row - 1)
        values(row:) = this%m_values(row + 1:)
        call move_alloc(values, this%m_values)
    end subroutine

! ******************************************************************************
! STRING_COLUMN MEMBERS
! ------------------------------------------------------------------------------
    pure function sc_init_empty(n) result(rst)
        !! Constructs an empty [[integer_column]] of the requested size.
        integer(int32), intent(in) :: n
            !! The number of rows to create.
        type(string_column) :: rst
            !! The newly constructed [[string_column]] object.
        
        if (n < 0) error stop FS_INVALID_INPUT_ERROR
        allocate(rst%m_values(n))
    end function

! ------------------------------------------------------------------------------
    pure function sc_init_array(x) result(rst)
        !! Constructs a new [[string_column]] by copying the contents of the
        !! supplied array.
        type(string), intent(in), dimension(:) :: x
            !! The array to copy.
        type(string_column) :: rst
            !! The newly constructed [[string_column]] object.

        allocate(rst%m_values(size(x)), source = x)
    end function

! ------------------------------------------------------------------------------
    pure function sc_count(this) result(rst)
        !! Gets the number of values stored in this column.
        class(string_column), intent(in) :: this
            !! The [[string_column]] object.
        integer(int32) :: rst
            !! The number of values.

        if (allocated(this%m_values)) then
            rst = size(this%m_values)
        else
            rst = 0
        end if
    end function

! ------------------------------------------------------------------------------
    pure function sc_get_type(this) result(rst)
        !! Gets the type stored by this column.
        class(string_column), intent(in) :: this
            !! The [[string_column]] object.
        integer(int32) :: rst
            !! The type stored by this column.
        rst = FS_COLUMN_STRING
    end function

! ------------------------------------------------------------------------------
    pure function sc_get(this, row) result(rst)
        !! Gets the requested value from the column.
        class(string_column), intent(in) :: this
            !! The [[string_column]] object.
        integer(int32), intent(in) :: row
            !! The row index of the value to retrieve.
        character(len = :), allocatable :: rst
            !! The value stored at the requested row.

        if (.not.allocated(this%m_values)) error stop FS_UNINITIALIZED_OBJECT_ERROR
        if (row < 1 .or. row > this%count()) error stop FS_INVALID_INPUT_ERROR
        rst = char(this%m_values(row))
    end function

! --------------------
    subroutine sc_set(this, row, x)
        !! Sets the requested value into the column.
        class(string_column), intent(inout) :: this
            !! The [[string_column]] object.
        integer(int32), intent(in) :: row
            !! The row index.
        character(len = *), intent(in) :: x
            !! The value stored at the requested row.

        if (.not.allocated(this%m_values)) error stop FS_UNINITIALIZED_OBJECT_ERROR
        if (row < 1 .or. row > this%count()) error stop FS_INVALID_INPUT_ERROR
        this%m_values(row) = x
    end subroutine

! ------------------------------------------------------------------------------
    pure subroutine sc_copy(this, dest)
        !! Copies the contents of this column to the supplied destination array.
        class(string_column), intent(in) :: this
            !! The [[string_column]] object.
        class(*), intent(inout), dimension(:) :: dest
            !! The destination array.

        if (size(dest) /= this%count()) error stop FS_ARRAY_SIZE_ERROR
        select type (dest)
        type is (string)
            if (allocated(this%m_values)) dest = this%m_values
        class default
            error stop FS_INVALID_ARGUMENT_ERROR
        end select
    end subroutine

! ------------------------------------------------------------------------------
    function sc_view(this) result(rst)
        !! Returns a borrowed, mutable view of this column's array.
        class(string_column), target, intent(inout) :: this
            !! The [[string_column]] object.
        class(*), pointer, dimension(:) :: rst
            !! The array pointer, or null if unallocated.

        nullify(rst)
        if (allocated(this%m_values)) rst => this%m_values
    end function

! ------------------------------------------------------------------------------
    function sc_slice(this, first, last) result(rst)
        !! Returns a borrowed, mutable slice of this column's array.
        class(string_column), target, intent(inout) :: this
            !! The [[string_column]] object.
        integer(int32), intent(in) :: first
            !! The one-based starting index of the slice.
        integer(int32), intent(in) :: last
            !! The one-based ending index of the slice.
        class(*), pointer, dimension(:) :: rst
            !! The array pointer, or null if unallocated.

        nullify(rst)
        if (first > last) error stop FS_INVALID_ARGUMENT_ERROR
        if (first < 1) error stop FS_INVALID_INPUT_ERROR
        if (last > this%count()) error stop FS_INVALID_INPUT_ERROR
        if (allocated(this%m_values)) rst => this%m_values(first:last)
    end function

! ------------------------------------------------------------------------------
    function sc_to_string(this, row, fmt) result(rst)
        !! Returns the value in the requested row as a string.
        class(string_column), intent(in) :: this
            !! The [[string_column]] object.
        integer(int32), intent(in) :: row
            !! The row index.
        character(len = *), intent(in), optional :: fmt
            !! The formatting string.
        character(len = :), allocatable :: rst
            !! The resulting string.
        
        if (.not.allocated(this%m_values)) error stop FS_UNINITIALIZED_OBJECT_ERROR
        if (row < 1 .or. row > this%count()) error stop FS_INVALID_INPUT_ERROR

        rst = char(this%m_values(row))
    end function

! ------------------------------------------------------------------------------
    subroutine sc_append(this, x)
        class(string_column), intent(inout) :: this
        class(*), intent(in) :: x

        call this%insert(this%count() + 1, x)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine sc_insert(this, row, x)
        class(string_column), intent(inout) :: this
        integer(int32), intent(in) :: row
        class(*), intent(in) :: x
        type(string), allocatable, dimension(:) :: values
        type(string) :: value
        integer(int32) :: n

        n = this%count()
        if (row < 1 .or. row > n + 1) error stop FS_INVALID_INPUT_ERROR
        select type (x)
        type is (string)
            value = x
        type is (character(len = *))
            value = x
        class default
            error stop FS_INVALID_ARGUMENT_ERROR
        end select
        allocate(values(n + 1))
        if (allocated(this%m_values)) then
            values(:row - 1) = this%m_values(:row - 1)
            values(row + 1:) = this%m_values(row:)
        end if
        values(row) = value
        call move_alloc(values, this%m_values)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine sc_remove(this, row)
        class(string_column), intent(inout) :: this
        integer(int32), intent(in) :: row
        type(string), allocatable, dimension(:) :: values
        integer(int32) :: n

        n = this%count()
        if (row < 1 .or. row > n) error stop FS_INVALID_INPUT_ERROR
        allocate(values(n - 1))
        values(:row - 1) = this%m_values(:row - 1)
        values(row:) = this%m_values(row + 1:)
        call move_alloc(values, this%m_values)
    end subroutine

! ******************************************************************************
! ADDITIONAL ROUTINES
! ------------------------------------------------------------------------------
    pure function dc_size(col) result(rst)
        !! Returns the number of items in the column.
        class(data_column), intent(in) :: col
            !! The [[data_column]] object.
        integer :: rst
            !! The size of the column.
        rst = col%count()
    end function

! ------------------------------------------------------------------------------
    function trim_trailing_zeros(value) result(rst)
        !! Trims trailing zeros from a string.
        character(len = *), intent(in) :: value
            !! The string on which to operate.
        character(len = :), allocatable :: rst
            !! The trimmed string.

        integer(int32) :: last_character, exponent_start

        rst = trim(value)
        exponent_start = scan(rst, "EeDd")
        if (exponent_start == 0) exponent_start = len(rst) + 1
        last_character = exponent_start - 1
        if (index(rst(:last_character), ".") == 0) return
        do while (last_character > 0)
            if (rst(last_character:last_character) /= "0") exit
            last_character = last_character - 1
        end do
        if (last_character > 0) then
            if (rst(last_character:last_character) == ".") last_character = last_character - 1
        end if
        rst = rst(:last_character) // rst(exponent_start:)
    end function

! ******************************************************************************
! REAL_COLUMN MATH OPERATIONS
! ------------------------------------------------------------------------------
    ! Element-wise math and arithmetic return independent, unnamed columns.
    ! Unallocated inputs are treated as empty for math and reductions.
    ! Empty reductions follow intrinsic identities: sum = 0, product = 1,
    ! maxval = -huge, and minval = huge. Domain errors follow intrinsic behavior.
    ! IEEE queries require all elements to match; allocated empty arrays return
    ! true, while unallocated columns return false.

    pure function rc_ieee_is_finite(x) result(rst)
        !! Overloads `ieee_is_finite` to test whether every column value is finite.
        !! An allocated empty column returns true; unallocated storage returns false.
        class(real_column), intent(in) :: x
            !! The [[real_column]] object.
        logical :: rst
            !! True if every value in the column is finite; else, false.

        if (allocated(x%m_values)) then
            rst = all(ieee_is_finite(x%m_values))
        else
            rst = .false.
        end if
    end function

! ------------------------------------------------------------------------------
    pure function rc_abs(x) result(rst)
        !! Overloads `abs` to compute the absolute value of each column element.
        class(real_column), intent(in) :: x
            !! The input [[real_column]]; unallocated storage is treated as empty.
        type(real_column) :: rst
            !! An independent, unnamed column with the same count as x.
        rst = real_column(0)
        if (allocated(x%m_values)) rst%m_values = abs(x%m_values)
    end function
! ------------------------------------------------------------------------------
    pure function rc_sin(x) result(rst)
        !! Overloads `sin` to compute the sine of each column element.
        class(real_column), intent(in) :: x
            !! The input angles in radians; unallocated storage is treated as empty.
        type(real_column) :: rst
            !! An independent, unnamed column of sine values with the same count as x.
        rst = real_column(0)
        if (allocated(x%m_values)) rst%m_values = sin(x%m_values)
    end function

! ------------------------------------------------------------------------------
    pure function rc_cos(x) result(rst)
        !! Overloads `cos` to compute the cosine of each column element.
        class(real_column), intent(in) :: x
            !! The input angles in radians; unallocated storage is treated as empty.
        type(real_column) :: rst
            !! An independent, unnamed column of cosine values with the same count as x.
        rst = real_column(0)
        if (allocated(x%m_values)) rst%m_values = cos(x%m_values)
    end function

! ------------------------------------------------------------------------------
    pure function rc_tan(x) result(rst)
        !! Overloads `tan` to compute the tangent of each column element.
        !! Singularities and nonfinite inputs follow intrinsic `tan` behavior.
        class(real_column), intent(in) :: x
            !! The input angles in radians; unallocated storage is treated as empty.
        type(real_column) :: rst
            !! An independent, unnamed column of tangent values with the same count as x.
        rst = real_column(0)
        if (allocated(x%m_values)) rst%m_values = tan(x%m_values)
    end function

! ------------------------------------------------------------------------------
    pure function rc_asin(x) result(rst)
        !! Overloads `asin` to compute the inverse sine of each column element.
        !! Values outside [-1, 1] follow intrinsic domain-error behavior.
        class(real_column), intent(in) :: x
            !! The input [[real_column]]; unallocated storage is treated as empty.
        type(real_column) :: rst
            !! An independent, unnamed column of angles in radians, with x's count.
        rst = real_column(0)
        if (allocated(x%m_values)) rst%m_values = asin(x%m_values)
    end function

! ------------------------------------------------------------------------------
    pure function rc_acos(x) result(rst)
        !! Overloads `acos` to compute the inverse cosine of each column element.
        !! Values outside [-1, 1] follow intrinsic domain-error behavior.
        class(real_column), intent(in) :: x
            !! The input [[real_column]]; unallocated storage is treated as empty.
        type(real_column) :: rst
            !! An independent, unnamed column of angles in radians, with x's count.
        rst = real_column(0)
        if (allocated(x%m_values)) rst%m_values = acos(x%m_values)
    end function

! ------------------------------------------------------------------------------
    pure function rc_atan(x) result(rst)
        !! Overloads the one-argument `atan` to compute each element's inverse tangent.
        class(real_column), intent(in) :: x
            !! The input [[real_column]]; unallocated storage is treated as empty.
        type(real_column) :: rst
            !! An independent, unnamed column of angles in radians, with x's count.
        rst = real_column(0)
        if (allocated(x%m_values)) rst%m_values = atan(x%m_values)
    end function

! ------------------------------------------------------------------------------
    pure function rc_exp(x) result(rst)
        !! Overloads `exp` to compute the natural exponential of each column element.
        !! Overflow and underflow follow intrinsic `exp` behavior.
        class(real_column), intent(in) :: x
            !! The input [[real_column]]; unallocated storage is treated as empty.
        type(real_column) :: rst
            !! An independent, unnamed column of exponential values, with x's count.
        rst = real_column(0)
        if (allocated(x%m_values)) rst%m_values = exp(x%m_values)
    end function

! ------------------------------------------------------------------------------
    pure function rc_log(x) result(rst)
        !! Overloads `log` to compute the natural logarithm of each column element.
        !! Nonpositive values follow intrinsic domain-error behavior.
        class(real_column), intent(in) :: x
            !! The input [[real_column]]; unallocated storage is treated as empty.
        type(real_column) :: rst
            !! An independent, unnamed column of logarithms with the same count as x.
        rst = real_column(0)
        if (allocated(x%m_values)) rst%m_values = log(x%m_values)
    end function

! ------------------------------------------------------------------------------
    pure function rc_log10(x) result(rst)
        !! Overloads `log10` to compute the base-10 logarithm of each column element.
        !! Nonpositive values follow intrinsic domain-error behavior.
        class(real_column), intent(in) :: x
            !! The input [[real_column]]; unallocated storage is treated as empty.
        type(real_column) :: rst
            !! An independent, unnamed column of logarithms with the same count as x.
        rst = real_column(0)
        if (allocated(x%m_values)) rst%m_values = log10(x%m_values)
    end function

! ------------------------------------------------------------------------------
    pure function rc_sqrt(x) result(rst)
        !! Overloads `sqrt` to compute the square root of each column element.
        !! Negative values follow intrinsic domain-error behavior.
        class(real_column), intent(in) :: x
            !! The input [[real_column]]; unallocated storage is treated as empty.
        type(real_column) :: rst
            !! An independent, unnamed column of square roots with x's count.
        rst = real_column(0)
        if (allocated(x%m_values)) rst%m_values = sqrt(x%m_values)
    end function

! ------------------------------------------------------------------------------
    pure function rc_sum(x) result(rst)
        !! Overloads `sum` to add all column values, without copying the array.
        class(real_column), intent(in) :: x
            !! The input [[real_column]]; unallocated storage is treated as empty.
        real(real64) :: rst
            !! The sum, or zero for an empty column. Uses intrinsic summation semantics.
        rst = 0.0d0
        if (allocated(x%m_values)) rst = sum(x%m_values)
    end function

! ------------------------------------------------------------------------------
    pure function rc_product(x) result(rst)
        !! Overloads `product` to multiply all column values, without copying the array.
        class(real_column), intent(in) :: x
            !! The input [[real_column]]; unallocated storage is treated as empty.
        real(real64) :: rst
            !! The product, or one for an empty column. Uses intrinsic product semantics.
        rst = 1.0d0
        if (allocated(x%m_values)) rst = product(x%m_values)
    end function

! ------------------------------------------------------------------------------
    pure function rc_maxval(x) result(rst)
        !! Overloads `maxval` to find the largest column value without copying the array.
        !! Nonfinite values follow intrinsic `maxval` behavior.
        class(real_column), intent(in) :: x
            !! The input [[real_column]]; unallocated storage is treated as empty.
        real(real64) :: rst
            !! The maximum, or -huge(0.0d0) for an empty column.
        rst = -huge(0.0d0)
        if (allocated(x%m_values)) rst = maxval(x%m_values)
    end function

! ------------------------------------------------------------------------------
    pure function rc_minval(x) result(rst)
        !! Overloads `minval` to find the smallest column value without copying the array.
        !! Nonfinite values follow intrinsic `minval` behavior.
        class(real_column), intent(in) :: x
            !! The input [[real_column]]; unallocated storage is treated as empty.
        real(real64) :: rst
            !! The minimum, or huge(0.0d0) for an empty column.
        rst = huge(0.0d0)
        if (allocated(x%m_values)) rst = minval(x%m_values)
    end function

! ------------------------------------------------------------------------------
    pure function rc_ieee_is_nan(x) result(rst)
        !! Overloads `ieee_is_nan` to test whether every column value is NaN.
        !! This is an all-elements query, not a test for the presence of any NaN.
        class(real_column), intent(in) :: x
            !! The [[real_column]] to inspect, without copying its array.
        logical :: rst
            !! True if all values are NaN, including an allocated empty array;
            !! false for unallocated storage.
        rst = .false.
        if (allocated(x%m_values)) rst = all(ieee_is_nan(x%m_values))
    end function

! ------------------------------------------------------------------------------
    pure function rc_ieee_is_negative(x) result(rst)
        !! Overloads `ieee_is_negative` to test the sign of every column value.
        !! Uses IEEE sign classification, including negative zero, rather than x < 0.
        class(real_column), intent(in) :: x
            !! The [[real_column]] to inspect, without copying its array.
        logical :: rst
            !! True if all values satisfy intrinsic `ieee_is_negative`, including
            !! an allocated empty array; false for unallocated storage.
        rst = .false.
        if (allocated(x%m_values)) rst = all(ieee_is_negative(x%m_values))
    end function

! ------------------------------------------------------------------------------
    pure function rc_ieee_is_normal(x) result(rst)
        !! Overloads `ieee_is_normal` to classify every column value.
        !! Zero, subnormal, infinite, and NaN values follow intrinsic classification.
        class(real_column), intent(in) :: x
            !! The [[real_column]] to inspect, without copying its array.
        logical :: rst
            !! True if all values satisfy intrinsic `ieee_is_normal`, including
            !! an allocated empty array; false for unallocated storage.
        rst = .false.
        if (allocated(x%m_values)) rst = all(ieee_is_normal(x%m_values))
    end function

! ------------------------------------------------------------------------------
    pure function rc_positive(x) result(rst)
        !! Overloads unary `+` to create a copy of the column's values.
        class(real_column), intent(in) :: x
            !! The operand [[real_column]]; unallocated storage is treated as empty.
        type(real_column) :: rst
            !! An independent, unnamed column with the same values and count as x.
        rst = real_column(0)
        if (allocated(x%m_values)) rst%m_values = +x%m_values
    end function

! ------------------------------------------------------------------------------
    pure function rc_negative(x) result(rst)
        !! Overloads unary `-` to negate each column value.
        class(real_column), intent(in) :: x
            !! The operand [[real_column]]; unallocated storage is treated as empty.
        type(real_column) :: rst
            !! An independent, unnamed column with negated values and x's count.
        rst = real_column(0)
        if (allocated(x%m_values)) rst%m_values = -x%m_values
    end function

! ------------------------------------------------------------------------------
    pure function rc_add(x, y) result(rst)
        !! Overloads `+` to add corresponding elements of two columns.
        !! Unequal counts cause ERROR STOP with [[FS_ARRAY_SIZE_ERROR]].
        class(real_column), intent(in) :: x, y
            !! The left and right operands; unallocated storage is treated as empty.
        type(real_column) :: rst
            !! An independent, unnamed column of x + y, with the operands' count.
        if (x%count() /= y%count()) error stop FS_ARRAY_SIZE_ERROR
        rst = real_column(0)
        if (x%count() > 0) rst%m_values = x%m_values + y%m_values
    end function

! ------------------------------------------------------------------------------
    pure function rc_subtract(x, y) result(rst)
        !! Overloads `-` to subtract corresponding elements of two columns.
        !! Unequal counts cause ERROR STOP with [[FS_ARRAY_SIZE_ERROR]].
        class(real_column), intent(in) :: x, y
            !! The left and right operands; unallocated storage is treated as empty.
        type(real_column) :: rst
            !! An independent, unnamed column of x - y, with the operands' count.
        if (x%count() /= y%count()) error stop FS_ARRAY_SIZE_ERROR
        rst = real_column(0)
        if (x%count() > 0) rst%m_values = x%m_values - y%m_values
    end function

! ------------------------------------------------------------------------------
    pure function rc_multiply(x, y) result(rst)
        !! Overloads `*` to multiply corresponding elements, not compute a dot product.
        !! Unequal counts cause ERROR STOP with [[FS_ARRAY_SIZE_ERROR]].
        class(real_column), intent(in) :: x, y
            !! The left and right operands; unallocated storage is treated as empty.
        type(real_column) :: rst
            !! An independent, unnamed column of x * y, with the operands' count.
        if (x%count() /= y%count()) error stop FS_ARRAY_SIZE_ERROR
        rst = real_column(0)
        if (x%count() > 0) rst%m_values = x%m_values * y%m_values
    end function

! ------------------------------------------------------------------------------
    pure function rc_divide(x, y) result(rst)
        !! Overloads `/` to divide corresponding elements of two columns.
        !! Unequal counts cause ERROR STOP with [[FS_ARRAY_SIZE_ERROR]].
        !! Zero divisors and nonfinite values follow intrinsic division behavior.
        class(real_column), intent(in) :: x, y
            !! The numerator and denominator columns; unallocated storage is empty.
        type(real_column) :: rst
            !! An independent, unnamed column of x / y, with the operands' count.
        if (x%count() /= y%count()) error stop FS_ARRAY_SIZE_ERROR
        rst = real_column(0)
        if (x%count() > 0) rst%m_values = x%m_values / y%m_values
    end function

! ------------------------------------------------------------------------------
    pure function rc_add_scalar(x, y) result(rst)
        !! Overloads `+` to add a scalar to each column value.
        class(real_column), intent(in) :: x
            !! The left operand [[real_column]]; unallocated storage is empty.
        real(real64), intent(in) :: y
            !! The scalar right operand, broadcast to every row.
        type(real_column) :: rst
            !! An independent, unnamed column of x + y, with x's count.
        rst = real_column(0)
        if (allocated(x%m_values)) rst%m_values = x%m_values + y
    end function

! ------------------------------------------------------------------------------
    pure function scalar_add_rc(x, y) result(rst)
        !! Overloads `+` to add each column value to a scalar.
        real(real64), intent(in) :: x
            !! The scalar left operand, broadcast to every row.
        class(real_column), intent(in) :: y
            !! The right operand [[real_column]]; unallocated storage is empty.
        type(real_column) :: rst
            !! An independent, unnamed column of x + y, with y's count.
        rst = y + x
    end function

! ------------------------------------------------------------------------------
    pure function rc_subtract_scalar(x, y) result(rst)
        !! Overloads `-` to subtract a scalar from each column value.
        class(real_column), intent(in) :: x
            !! The left operand [[real_column]]; unallocated storage is empty.
        real(real64), intent(in) :: y
            !! The scalar right operand, broadcast to every row.
        type(real_column) :: rst
            !! An independent, unnamed column of x - y, with x's count.
        rst = real_column(0)
        if (allocated(x%m_values)) rst%m_values = x%m_values - y
    end function

! ------------------------------------------------------------------------------
    pure function scalar_subtract_rc(x, y) result(rst)
        !! Overloads `-` to subtract each column value from a scalar.
        real(real64), intent(in) :: x
            !! The scalar left operand, broadcast to every row.
        class(real_column), intent(in) :: y
            !! The right operand [[real_column]]; unallocated storage is empty.
        type(real_column) :: rst
            !! An independent, unnamed column of x - y, with y's count.
        rst = real_column(0)
        if (allocated(y%m_values)) rst%m_values = x - y%m_values
    end function

! ------------------------------------------------------------------------------
    pure function rc_multiply_scalar(x, y) result(rst)
        !! Overloads `*` to multiply each column value by a scalar.
        class(real_column), intent(in) :: x
            !! The left operand [[real_column]]; unallocated storage is empty.
        real(real64), intent(in) :: y
            !! The scalar right operand, broadcast to every row.
        type(real_column) :: rst
            !! An independent, unnamed column of x * y, with x's count.
        rst = real_column(0)
        if (allocated(x%m_values)) rst%m_values = x%m_values * y
    end function

! ------------------------------------------------------------------------------
    pure function scalar_multiply_rc(x, y) result(rst)
        !! Overloads `*` to multiply a scalar by each column value.
        real(real64), intent(in) :: x
            !! The scalar left operand, broadcast to every row.
        class(real_column), intent(in) :: y
            !! The right operand [[real_column]]; unallocated storage is empty.
        type(real_column) :: rst
            !! An independent, unnamed column of x * y, with y's count.
        rst = y * x
    end function

! ------------------------------------------------------------------------------
    pure function rc_divide_scalar(x, y) result(rst)
        !! Overloads `/` to divide each column value by a scalar.
        !! Zero divisors and nonfinite values follow intrinsic division behavior.
        class(real_column), intent(in) :: x
            !! The numerator [[real_column]]; unallocated storage is empty.
        real(real64), intent(in) :: y
            !! The scalar denominator, broadcast to every row.
        type(real_column) :: rst
            !! An independent, unnamed column of x / y, with x's count.
        rst = real_column(0)
        if (allocated(x%m_values)) rst%m_values = x%m_values / y
    end function

! ------------------------------------------------------------------------------
    pure function scalar_divide_rc(x, y) result(rst)
        !! Overloads `/` to divide a scalar by each column value.
        !! Zero divisors and nonfinite values follow intrinsic division behavior.
        real(real64), intent(in) :: x
            !! The scalar numerator, broadcast to every row.
        class(real_column), intent(in) :: y
            !! The denominator [[real_column]]; unallocated storage is empty.
        type(real_column) :: rst
            !! An independent, unnamed column of x / y, with y's count.
        rst = real_column(0)
        if (allocated(y%m_values)) rst%m_values = x / y%m_values
    end function

! ------------------------------------------------------------------------------
end module
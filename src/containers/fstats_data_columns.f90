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
        procedure, public :: to_string => sc_to_string
        procedure, public :: append => sc_append
        procedure, public :: insert => sc_insert
        procedure, public :: remove => sc_remove
    end type

    interface string_column
        module procedure :: sc_init_empty
        module procedure :: sc_init_array
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

! ******************************************************************************
! ADDITIONAL ROUTINES
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

! ------------------------------------------------------------------------------
end module
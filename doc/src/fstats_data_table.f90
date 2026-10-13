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

module fstats_data_table
    use iso_fortran_env, only : int32
    use fstats_data_columns
    use fstats_errors
    implicit none
    private
    public :: data_table
    public :: data_table_view
    public :: assignment(=)

    type :: data_table_column_slot
        !! Owns one individually allocated polymorphic data column.
        class(data_column), pointer :: value => null()
            !! The column object stored in this slot, or null when unused.
    end type

    type :: data_table_column_collection
        !! Dynamically sized owner for the table's heterogeneous columns.
        type(data_table_column_slot), allocatable, dimension(:) :: slots
            !! Storage for column handles; capacity can exceed the active count.
        integer(int32) :: active_count = 0
            !! Number of initialized column slots.
    contains
        procedure, private :: reserve => dtcc_reserve
        procedure, public :: count => dtcc_count
        procedure, public :: append => dtcc_append
        procedure, public :: insert => dtcc_insert
        procedure, public :: remove => dtcc_remove
        procedure, public :: get => dtcc_get
        procedure, public :: clear => dtcc_clear
        procedure, public :: take_contents => dtcc_take_contents
        final :: dtcc_finalize
    end type

    type :: data_table_state
        !! Shared mutable state retained by its table and all live views.
        type(data_table_column_collection) :: columns
            !! The owned heterogeneous columns.
        integer(int32) :: reference_count = 1
            !! Number of table and view handles retaining this state.
    end type

    type data_table
        !! Defines a table of data.
        type(data_table_state), pointer, private :: m_state => null()
            !! Shared mutable state; null represents a table with no allocation.
    contains
        procedure, public :: row_count => dt_row_count
        procedure, public :: column_count => dt_column_count
        procedure, public :: add_column => dt_add_column
        procedure, public :: remove_column => dt_remove_column
        procedure, public :: insert_column => dt_insert_column
        procedure, public :: get_column => dt_get_column
        procedure, public :: view => dt_view
        procedure, public :: slice => dt_slice
        final :: dt_finalize
    end type

    type data_table_view
        !! A live positional selection of rows and columns in shared table state.
        type(data_table_state), pointer, private :: m_state => null()
            !! Shared mutable state retained independently of the source table.
        integer(int32), private :: row_first = 1
            !! First selected one-based row position.
        integer(int32), private :: row_last = 0
            !! Last selected row; zero selects through the current final row.
        integer(int32), private :: column_first = 1
            !! First selected one-based column position.
        integer(int32), private :: column_last = 0
            !! Last selected column; zero selects through the current final column.
    contains
        procedure, public :: row_count => dtv_row_count
        procedure, public :: column_count => dtv_column_count
        procedure, public :: get_column => dtv_get_column
        procedure, public :: column_slice => dtv_column_slice
        final :: dtv_finalize
    end type

    interface assignment(=)
        !! Defines deep-copy table assignment and shared view assignment.
        module procedure :: dt_assign
        module procedure :: dtv_assign
    end interface

contains
! ------------------------------------------------------------------------------
    pure function dtcc_count(this) result(rst)
        !! Returns the number of active column slots.
        class(data_table_column_collection), intent(in) :: this
            !! The collection being queried.
        integer(int32) :: rst
            !! The number of stored columns.
        rst = this%active_count
    end function

! ------------------------------------------------------------------------------
    subroutine dtcc_reserve(this, needed)
        !! Grows column-handle storage without moving the separately allocated columns.
        class(data_table_column_collection), intent(inout) :: this
            !! The collection whose capacity may need to grow.
        integer(int32), intent(in) :: needed
            !! The minimum required number of slots.

        type(data_table_column_slot), allocatable, dimension(:) :: expanded
        integer(int32) :: capacity, new_capacity, i

        capacity = 0
        if (allocated(this%slots)) capacity = size(this%slots)
        if (capacity >= needed) return

        new_capacity = max(4, capacity * 2, needed)
        allocate(expanded(new_capacity))
        do i = 1, this%active_count
            expanded(i)%value => this%slots(i)%value
            nullify(this%slots(i)%value)
        end do
        call move_alloc(expanded, this%slots)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine dtcc_append(this, column)
        !! Copies a column into the next collection slot.
        class(data_table_column_collection), intent(inout) :: this
            !! The collection receiving the column.
        class(data_column), intent(in) :: column
            !! The source column to copy.

        class(data_column), pointer :: owned_column

        allocate(owned_column, source = column)
        call this%reserve(this%active_count + 1)
        this%active_count = this%active_count + 1
        this%slots(this%active_count)%value => owned_column
    end subroutine

! ------------------------------------------------------------------------------
    subroutine dtcc_insert(this, index, column)
        !! Copies a column into the requested position while preserving other column objects.
        class(data_table_column_collection), intent(inout) :: this
            !! The collection receiving the column.
        integer(int32), intent(in) :: index
            !! One-based insertion position, from one through count plus one.
        class(data_column), intent(in) :: column
            !! The source column to copy.

        class(data_column), pointer :: owned_column
        integer(int32) :: i

        if (index < 1 .or. index > this%active_count + 1) then
            error stop FS_INVALID_INPUT_ERROR
        end if
        allocate(owned_column, source = column)
        call this%reserve(this%active_count + 1)
        do i = this%active_count, index, -1
            this%slots(i + 1)%value => this%slots(i)%value
            nullify(this%slots(i)%value)
        end do
        this%slots(index)%value => owned_column
        this%active_count = this%active_count + 1
    end subroutine

! ------------------------------------------------------------------------------
    subroutine dtcc_remove(this, index)
        !! Removes and deallocates a column, shifting remaining handles in place.
        class(data_table_column_collection), intent(inout) :: this
            !! The collection from which the column is removed.
        integer(int32), intent(in) :: index
            !! One-based position of the column to remove.

        integer(int32) :: i

        if (this%active_count == 0) return
        if (index < 1 .or. index > this%active_count) then
            error stop FS_INVALID_INPUT_ERROR
        end if
        deallocate(this%slots(index)%value)
        nullify(this%slots(index)%value)
        do i = index, this%active_count - 1
            this%slots(i)%value => this%slots(i + 1)%value
            nullify(this%slots(i + 1)%value)
        end do
        this%active_count = this%active_count - 1
    end subroutine

! ------------------------------------------------------------------------------
    subroutine dtcc_clear(this)
        !! Deallocates all active columns while retaining collection capacity.
        class(data_table_column_collection), intent(inout) :: this
            !! The collection to empty.

        integer(int32) :: i

        do i = 1, this%active_count
            if (associated(this%slots(i)%value)) then
                deallocate(this%slots(i)%value)
                nullify(this%slots(i)%value)
            end if
        end do
        this%active_count = 0
    end subroutine

! ------------------------------------------------------------------------------
    subroutine dtcc_take_contents(this, source)
        !! Replaces this collection's contents by transferring owned column handles.
        class(data_table_column_collection), intent(inout) :: this
            !! The collection receiving the transferred contents.
        class(data_table_column_collection), intent(inout) :: source
            !! The collection that relinquishes its contents.

        call this%clear()
        call move_alloc(source%slots, this%slots)
        this%active_count = source%active_count
        source%active_count = 0
    end subroutine

! ------------------------------------------------------------------------------
    function dtcc_get(this, index) result(rst)
        !! Returns a borrowed pointer to a column in the collection.
        class(data_table_column_collection), intent(in) :: this
            !! The collection being queried.
        integer(int32), intent(in) :: index
            !! One-based position of the requested column.
        class(data_column), pointer :: rst
            !! The stored column; invalid indices cause an error stop.

        nullify(rst)
        if (index < 1 .or. index > this%active_count) then
            error stop FS_INVALID_INPUT_ERROR
        end if
        rst => this%slots(index)%value
    end function

! ------------------------------------------------------------------------------
    subroutine dtcc_finalize(this)
        !! Deallocates every column still owned by the collection.
        type(data_table_column_collection), intent(inout) :: this
            !! The collection being finalized.

        integer(int32) :: i

        do i = 1, this%active_count
            if (associated(this%slots(i)%value)) then
                deallocate(this%slots(i)%value)
                nullify(this%slots(i)%value)
            end if
        end do
        this%active_count = 0
    end subroutine

! ------------------------------------------------------------------------------
    subroutine dt_ensure_state(this)
        !! Allocates initial shared state when a table is first used.
        class(data_table), intent(inout) :: this
            !! The table whose state is initialized.

        if (associated(this%m_state)) return
        allocate(this%m_state)
        this%m_state%reference_count = 1
    end subroutine

! ------------------------------------------------------------------------------
    subroutine dt_release_state(state)
        !! Releases one table or view reference to shared state.
        type(data_table_state), pointer, intent(inout) :: state
            !! The state reference to release; null on return.

        if (.not.associated(state)) return
        state%reference_count = state%reference_count - 1
        if (state%reference_count == 0) deallocate(state)
        nullify(state)
    end subroutine

! ------------------------------------------------------------------------------
    function dt_clone_state(source) result(clone)
        !! Deep-copies a table state for value-semantic table assignment.
        type(data_table_state), pointer, intent(in) :: source
            !! The state to clone, or null for an unallocated empty table.
        type(data_table_state), pointer :: clone
            !! The independent copy, or null when source is null.

        class(data_column), pointer :: column
        integer(int32) :: i

        nullify(clone)
        if (.not.associated(source)) return
        allocate(clone)
        clone%reference_count = 1
        do i = 1, source%columns%count()
            column => source%columns%get(i)
            call clone%columns%append(column)
        end do
    end function

! ------------------------------------------------------------------------------
    subroutine dt_assign(lhs, rhs)
        !! Deep-copies table contents while preserving state identity for attached views.
        type(data_table), intent(inout) :: lhs
            !! The table receiving the copied contents.
        type(data_table), intent(in) :: rhs
            !! The table whose contents are copied.

        type(data_table_state), pointer :: copied_state

        copied_state => dt_clone_state(rhs%m_state)
        if (associated(lhs%m_state)) then
            if (associated(copied_state)) then
                call lhs%m_state%columns%take_contents(copied_state%columns)
                deallocate(copied_state)
            else
                call lhs%m_state%columns%clear()
            end if
        else
            lhs%m_state => copied_state
        end if
    end subroutine

! ------------------------------------------------------------------------------
    subroutine dtv_assign(lhs, rhs)
        !! Assigns a view by sharing its state and copying its positional bounds.
        type(data_table_view), intent(inout) :: lhs
            !! The view receiving the shared state and selector.
        type(data_table_view), intent(in) :: rhs
            !! The view whose state and selector are copied.

        type(data_table_state), pointer :: shared_state
        integer(int32) :: row_first, row_last, column_first, column_last

        shared_state => rhs%m_state
        row_first = rhs%row_first
        row_last = rhs%row_last
        column_first = rhs%column_first
        column_last = rhs%column_last
        if (associated(shared_state)) then
            shared_state%reference_count = shared_state%reference_count + 1
        end if
        call dt_release_state(lhs%m_state)
        lhs%m_state => shared_state
        lhs%row_first = row_first
        lhs%row_last = row_last
        lhs%column_first = column_first
        lhs%column_last = column_last
    end subroutine

! ------------------------------------------------------------------------------
    subroutine dt_finalize(this)
        !! Releases the state retained by a table when it leaves scope.
        type(data_table), intent(inout) :: this
            !! The table being finalized.

        call dt_release_state(this%m_state)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine dtv_finalize(this)
        !! Releases the state retained by a view when it leaves scope.
        type(data_table_view), intent(inout) :: this
            !! The view being finalized.

        call dt_release_state(this%m_state)
    end subroutine

! ------------------------------------------------------------------------------
    pure function dt_bounded_count(first, last, extent) result(rst)
        !! Computes the visible size of a positional range within a current extent.
        integer(int32), intent(in) :: first
            !! First one-based position in the selector.
        integer(int32), intent(in) :: last
            !! Last selector position; zero denotes the current end.
        integer(int32), intent(in) :: extent
            !! Current number of available positions.
        integer(int32) :: rst
            !! Number of selected positions currently available.

        integer(int32) :: upper

        rst = 0
        if (extent < 1 .or. first > extent) return
        upper = extent
        if (last > 0) upper = min(last, extent)
        if (upper >= first) rst = upper - first + 1
    end function

! ------------------------------------------------------------------------------
    function dt_make_view(this, row_first, row_last, column_first, column_last) result(rst)
        !! Creates a retained view state with validated positional selectors.
        class(data_table), intent(inout) :: this
            !! The source table whose state is shared.
        integer(int32), intent(in) :: row_first
            !! First selected row position.
        integer(int32), intent(in) :: row_last
            !! Last selected row position, or zero for the current final row.
        integer(int32), intent(in) :: column_first
            !! First selected column position.
        integer(int32), intent(in) :: column_last
            !! Last selected column position, or zero for the current final column.
        type(data_table_view) :: rst
            !! The live view retaining the source table state.

        call dt_ensure_state(this)
        if (row_first < 1 .or. row_last < 0 .or. &
            (row_last > 0 .and. row_last < row_first)) then
            error stop FS_INVALID_INPUT_ERROR
        end if
        if (column_first < 1 .or. column_last < 0 .or. &
            (column_last > 0 .and. column_last < column_first)) then
            error stop FS_INVALID_INPUT_ERROR
        end if

        this%m_state%reference_count = this%m_state%reference_count + 1
        rst%m_state => this%m_state
        rst%row_first = row_first
        rst%row_last = row_last
        rst%column_first = column_first
        rst%column_last = column_last
    end function

! ------------------------------------------------------------------------------
    function dt_view(this, row_first, row_last, column_first, column_last) result(rst)
        !! Creates a live view using optional one-based row and column bounds.
        !! Bounds are reapplied to current positions and clipped to the current
        !! table shape, so structural edits can change the selected values.
        !! Omitted final bounds expand as rows or columns are appended; explicit
        !! final bounds remain fixed.
        class(data_table), intent(inout) :: this
            !! The source table whose shared state is retained.
        integer(int32), intent(in), optional :: row_first
            !! First selected row; defaults to one.
        integer(int32), intent(in), optional :: row_last
            !! Inclusive final row; omitted means through the current final row.
        integer(int32), intent(in), optional :: column_first
            !! First selected column; defaults to one.
        integer(int32), intent(in), optional :: column_last
            !! Inclusive final column; omitted means through the current final column.
        type(data_table_view) :: rst
            !! The live positional view.

        integer(int32) :: first_row, last_row, first_column, last_column

        first_row = 1
        last_row = 0
        first_column = 1
        last_column = 0
        if (present(row_first)) first_row = row_first
        if (present(row_last)) last_row = row_last
        if (present(column_first)) first_column = column_first
        if (present(column_last)) last_column = column_last
        if (present(row_last)) then
            if (row_last < 1) error stop FS_INVALID_INPUT_ERROR
        end if
        if (present(column_last)) then
            if (column_last < 1) error stop FS_INVALID_INPUT_ERROR
        end if
        rst = dt_make_view(this, first_row, last_row, first_column, last_column)
    end function

! ------------------------------------------------------------------------------
    function dt_slice(this, first, last) result(rst)
        !! Creates a live row slice spanning all current and future columns.
        class(data_table), intent(inout) :: this
            !! The source table whose shared state is retained.
        integer(int32), intent(in) :: first
            !! First selected one-based row position.
        integer(int32), intent(in) :: last
            !! Inclusive final row position.
        type(data_table_view) :: rst
            !! A view restricted to the requested row positions.

        if (first < 1 .or. last < first) error stop FS_INVALID_INPUT_ERROR
        rst = dt_make_view(this, first, last, 1, 0)
    end function

! ------------------------------------------------------------------------------
    function dt_row_count(this) result(rst)
        !! Gets the number of rows in the data table.
        class(data_table), intent(in) :: this
            !! The table being queried.
        integer(int32) :: rst
            !! The number of rows, or zero when there are no columns.

        rst = 0
        if (.not.associated(this%m_state)) return
        if (this%m_state%columns%count() == 0) return
        rst = this%m_state%columns%slots(1)%value%count()
    end function

! ------------------------------------------------------------------------------
    pure function dt_column_count(this) result(rst)
        !! Gets the number of columns in the data table.
        class(data_table), intent(in) :: this
            !! The table being queried.
        integer(int32) :: rst
            !! The number of columns.

        rst = 0
        if (associated(this%m_state)) rst = this%m_state%columns%count()
    end function

! ------------------------------------------------------------------------------
    subroutine dt_add_column(this, col)
        !! Appends an owned copy of a column to the table.
        class(data_table), intent(inout) :: this
            !! The table receiving the column.
        class(data_column), intent(in) :: col
            !! The column to copy; its row count must match the table.

        if (this%column_count() > 0) then
            if (col%count() /= this%row_count()) error stop FS_ARRAY_SIZE_ERROR
        end if
        call dt_ensure_state(this)
        call this%m_state%columns%append(col)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine dt_remove_column(this, col)
        !! Removes a column from the table and shifts later positions left.
        class(data_table), intent(inout) :: this
            !! The table from which the column is removed.
        integer(int32), intent(in) :: col
            !! One-based position of the column to remove.

        if (.not.associated(this%m_state)) return
        call this%m_state%columns%remove(col)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine dt_insert_column(this, index, col)
        !! Inserts an owned copy of a column at the requested position.
        class(data_table), intent(inout) :: this
            !! The table receiving the column.
        integer(int32), intent(in) :: index
            !! One-based insertion position, from one through count plus one.
        class(data_column), intent(in) :: col
            !! The column to copy; its row count must match the table.

        if (this%column_count() > 0) then
            if (col%count() /= this%row_count()) error stop FS_ARRAY_SIZE_ERROR
        end if
        call dt_ensure_state(this)
        call this%m_state%columns%insert(index, col)
    end subroutine

! ------------------------------------------------------------------------------
    function dt_get_column(this, col) result(rst)
        !! Gets a borrowed pointer to a stored column without copying it.
        !! The pointer survives insertion, but removal of that column or release
        !! of the final table/view owner invalidates it.
        class(data_table), intent(in) :: this
            !! The table being queried.
        integer(int32), intent(in) :: col
            !! One-based position of the requested column.
        class(data_column), pointer :: rst
            !! The stored column; invalid indices cause an error stop.

        nullify(rst)
        if (.not.associated(this%m_state)) error stop FS_INVALID_INPUT_ERROR
        rst => this%m_state%columns%get(col)
    end function

! ------------------------------------------------------------------------------
    function dtv_row_count(this) result(rst)
        !! Returns the number of selected rows currently available.
        class(data_table_view), intent(in) :: this
            !! The view being queried.
        integer(int32) :: rst
            !! The visible row count after clipping to the current table shape.

        integer(int32) :: extent

        extent = 0
        if (associated(this%m_state)) then
            if (this%m_state%columns%count() > 0) then
                extent = this%m_state%columns%slots(1)%value%count()
            end if
        end if
        rst = dt_bounded_count(this%row_first, this%row_last, extent)
    end function

! ------------------------------------------------------------------------------
    pure function dtv_column_count(this) result(rst)
        !! Returns the number of selected columns currently available.
        class(data_table_view), intent(in) :: this
            !! The view being queried.
        integer(int32) :: rst
            !! The visible column count after clipping to the current table shape.

        integer(int32) :: extent

        extent = 0
        if (associated(this%m_state)) extent = this%m_state%columns%count()
        rst = dt_bounded_count(this%column_first, this%column_last, extent)
    end function

! ------------------------------------------------------------------------------
    function dtv_get_column(this, col) result(rst)
        !! Gets a borrowed pointer to a selected column's full underlying column.
        !! The pointer is not restricted to the selected rows; use `column_slice`
        !! for those rows. Removing the column invalidates the pointer.
        class(data_table_view), intent(in) :: this
            !! The view selecting the column.
        integer(int32), intent(in) :: col
            !! One-based position within the view's selected columns.
        class(data_column), pointer :: rst
            !! The stored column; invalid view positions cause an error stop.

        integer(int32) :: actual_column

        nullify(rst)
        if (col < 1 .or. col > this%column_count()) error stop FS_INVALID_INPUT_ERROR
        actual_column = this%column_first + col - 1
        rst => this%m_state%columns%get(actual_column)
    end function

! ------------------------------------------------------------------------------
    function dtv_column_slice(this, col) result(rst)
        !! Borrows the selected row range from one selected column.
        !! The pointer is invalidated if the column is removed or its array is
        !! reallocated; reacquire it after such edits.
        class(data_table_view), intent(in) :: this
            !! The view selecting rows and columns.
        integer(int32), intent(in) :: col
            !! One-based position of the column within the view.
        class(*), pointer, dimension(:) :: rst
            !! Mutable selected values; reacquire after a source-column reallocation.

        class(data_column), pointer :: selected_column
        integer(int32) :: actual_column, actual_last_row

        nullify(rst)
        if (col < 1 .or. col > this%column_count()) error stop FS_INVALID_INPUT_ERROR
        if (this%row_count() == 0) return
        actual_column = this%column_first + col - 1
        actual_last_row = this%row_first + this%row_count() - 1
        selected_column => this%m_state%columns%get(actual_column)
        rst => selected_column%slice(this%row_first, actual_last_row)
    end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
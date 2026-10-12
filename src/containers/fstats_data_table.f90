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
    use iso_fortran_env
    use collections
    use fstats_data_columns
    use fstats_errors
    implicit none
    private
    public :: data_table
    
    type data_table
        !! Defines a table of data.
        type(list), private :: m_columns
    contains
        procedure, public :: row_count => dt_row_count
        procedure, public :: column_count => dt_column_count
        procedure, public :: add_column => dt_add_column
        procedure, public :: remove_column => dt_remove_column
        procedure, public :: insert_column => dt_insert_column
        procedure, public :: get_column => dt_get_column
    end type

contains
! ------------------------------------------------------------------------------
    function dt_row_count(this) result(rst)
        !! Gets the number of rows in the data table.
        class(data_table), intent(in) :: this
            !! The [[data_table]] object.
        integer(int32) :: rst
            !! The row count.
        
        class(data_column), pointer :: col
        if (this%column_count() == 0) then
            rst = 0
        else
            col => this%get_column(1)
            rst = col%count()
        end if
    end function

! ------------------------------------------------------------------------------
    pure function dt_column_count(this) result(rst)
        !! Gets the number of columns in the data table.
        class(data_table), intent(in) :: this
            !! The [[data_table]] object.
        integer(int32) :: rst
            !! The column count.
        rst = this%m_columns%count()
    end function

! ------------------------------------------------------------------------------
    subroutine dt_add_column(this, col)
        !! Appends a column onto the [[data_table]].
        class(data_table), intent(inout) :: this
            !! The [[data_table]] object.
        class(data_column), intent(in) :: col
            !! The column to add.  The columns is copied internally, and the 
            !! data table takes over management of allocated memory.  An 
            !! [[FS_ARRAY_SIZE_ERROR]] is thrown if the number of items in the
            !! supplied column does not match the number of rows in the 
            !! data table.

        if (this%column_count() > 0) then
            if (col%count() /= this%row_count()) error stop FS_ARRAY_SIZE_ERROR
        end if

        call this%m_columns%push(col, manage = .true.)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine dt_remove_column(this, col)
        !! Removes the requested column from the [[data_table]].
        class(data_table), intent(inout) :: this
            !! The [[data_table]] object.
        integer(int32), intent(in) :: col
            !! The index of the column to remove.
        call this%m_columns%remove(col)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine dt_insert_column(this, index, col)
        !! Inserts a column into the [[data_table]].
        class(data_table), intent(inout) :: this
            !! The [[data_table]] object.
        integer(int32), intent(in) :: index
            !! The index where the column is to be inserted.
        class(data_column), intent(in) :: col
            !! The column to insert.  An [[FS_ARRAY_SIZE_ERROR]] is thrown if 
            !! the number of items in the supplied column does not match the 
            !! number of rows in the data table.

        if (this%column_count() > 0) then
            if (col%count() /= this%row_count()) error stop FS_ARRAY_SIZE_ERROR
        end if
        
        call this%m_columns%insert(index, col, manage = .true.)
    end subroutine

! ------------------------------------------------------------------------------
    function dt_get_column(this, col) result(rst)
        !! Gets a borrowed pointer to the stored column without copying it.
        !! Changes through this pointer modify the table's column. The caller
        !! must not deallocate it; removing the column or destroying the table
        !! invalidates the pointer. Invalid indices cause the list's index error.
        class(data_table), intent(in) :: this
            !! The [[data_table]] object.
        integer(int32), intent(in) :: col
            !! The one-based index of the column to retrieve.
        class(data_column), pointer :: rst
            !! The requested column, retaining its concrete dynamic type;
            !! null if the stored object is unassociated or is not a data column.

        class(*), pointer :: item

        nullify(rst)
        item => this%m_columns%get(col)
        if (.not.associated(item)) return
        select type (item)
        class is (data_column)
            rst => item
        end select
    end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
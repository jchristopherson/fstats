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

program example
    use iso_fortran_env
    use fstats
    use fplot_core
    use csv_module
    implicit none

    ! Local Variables
    real(real64) :: dt
    real(real64), allocatable, dimension(:,:) :: v
    type(data_table) :: tbl

    type(plot_2d) :: plt
    type(plot_data_2d) :: pd
    class(plot_axis), pointer :: xAxis, yAxis
    class(data_column), pointer :: t, x

    ! Read in the data
    call read_delimited_table("examples/data/noise_data.csv", tbl, header_row = 1)
    t => tbl%get_column(1)
    x => tbl%get_column(2)
    select type (t)
    class is (real_column)
        dt = t%get(2) - t%get(1)
    class default
        dt = 0.0d0
    end select
    print "(A, F0.4)", "dt = ", dt

    ! Compute the Allan variance
    v = allan_variance(x, dt)

    ! Plot the results
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    call xAxis%set_is_log_scaled(.true.)
    call yAxis%set_is_log_scaled(.true.)
    call yAxis%set_use_default_tic_label_format(.false.)
    call yAxis%set_tic_label_format("%0.0e")
    call xAxis%set_title("Averaging Time")
    call yAxis%set_title("Deviation")
    call pd%define_data(v(:,1), sqrt(v(:,2)))
    call plt%push(pd)
    call plt%draw()
end program
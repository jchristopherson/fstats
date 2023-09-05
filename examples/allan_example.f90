program example
    use iso_fortran_env
    use fstats
    use fplot_core
    use csv_module
    implicit none

    ! Local Variables
    real(real64) :: dt
    real(real64), allocatable, dimension(:) :: t, x
    real(real64), allocatable, dimension(:,:) :: v
    logical :: ok
    type(csv_file) :: file
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd
    class(plot_axis), pointer :: xAxis, yAxis

    ! Read in the data
    call file%read("examples/data/noise_data.csv", header_row = 1, status_ok = ok)
    if (.not.ok) then
        print "(A)", "Could not open the data file."
        stop
    end if
    call file%get(1, t, ok)
    if (.not.ok) then
        print "(A)", "Could not extract the time data."
        stop
    end if
    call file%get(2, x, ok)
    if (.not.ok) then
        print "(A)", "Could not extract the signal data."
        stop
    end if
    dt = t(2) - t(1)

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
    call yAxis%set_title("Variance")
    call pd%define_data(v(:,1), v(:,2))
    call plt%push(pd)
    call plt%draw()
end program
program example
    use iso_fortran_env
    use fstats
    implicit none

    ! Local Variables
    integer(int32) :: i, vars(3), tbl(24, 3)

    ! Define the number of design points for each of the 3 factors to study
    vars = [2, 4, 3]

    ! Determine the design table
    call full_factorial(vars, tbl)

    ! Display the table
    do i = 1, size(tbl, 1)
        print *, tbl(i,:)
    end do
end program
module fstats_types
    use iso_fortran_env
    implicit none

    type array_container
        !! Provides a container for a real-valued array.  A practical use of
        !! this construct is in the construction of jagged arrays.
        real(real64), allocatable, dimension(:) :: x
            !! The array.
    end type
end module
module external

    interface
        subroutine SGEMM(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, &
            c, ldc)
            use iso_fortran_env, only : int32, real32
            character, intent(in) :: transa, transb
            integer(int32), intent(in) :: m, n, k, lda, ldb, ldc
            real(real32), intent(in) :: alpha, a(lda,*), b(ldb,*), beta
            real(real32), intent(inout) :: c(ldc,*)
        end subroutine
    end interface
    
end module
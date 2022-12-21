! fstats.f90

module fstats
    use iso_fortran_env
    implicit none
    private
    public :: beta
    public :: regularized_beta
    public :: incomplete_beta
    public :: digamma
    public :: incomplete_gamma_upper
    public :: incomplete_gamma_lower


! ******************************************************************************
! SPECIAL FUNCTIONS
! ------------------------------------------------------------------------------
    !> Computes the beta function.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! real(real64) function beta(real(real64) a, real(real64) b)
    !! real(real32) function beta(real(real32) a, real(real32) b)
    !! @endcode
    !!
    !! @param[in] a The first argument of the function.
    !! @param[in] b The second argument of the function.
    !!
    !! @return The value of the beta function at @p a and @p b.
    !!
    !! @remarks The beta function is related to the gamma function
    !! by the following relationship \f$ \beta(a,b) = 
    !! \frac{\Gamma(a) \Gamma(b)}{\Gamma(a + b)} \f$.
    interface beta
        module procedure :: beta_real64
        module procedure :: beta_real32
    end interface

    !> Computes the regularized beta function.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! real(real64) function regularized_beta(real(real64) a, real(real64) b, real(real64) x)
    !! real(real32) function regularized_beta(real(real32) a, real(real32) b, real(real32) x)
    !! @endcode
    !!
    !! @param[in] a The first argument of the function.
    !! @param[in] b The second argument of the function.
    !! @param[in] x The upper limit of the integration.
    !!
    !! @return The value of the incomplete beta function.
    !!
    !! @remarks The regularized beta function is defined as the ratio between
    !! the incomplete beta function and the beta function: \f$ I_{x}(a,b) = 
    !! \frac{\beta(x;a,b)}{\beta(a,b)} \f$.
    interface regularized_beta
        module procedure :: regularized_beta_real64
        module procedure :: regularized_beta_real32
    end interface

    !> Computes the incomplete beta function.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! real(real64) function incomplete_beta(real(real64) a, real(real64) b, real(real64) x)
    !! real(real32) function incomplete_beta(real(real32) a, real(real32) b, real(real32) x)
    !! @endcode
    !!
    !! @param[in] a The first argument of the function.
    !! @param[in] b The second argument of the function.
    !! @param[in] x The upper limit of the integration.
    !!
    !! @return The value of the incomplete beta function.
    !!
    !! @remarks The incomplete beta function is defind as \f$ \beta(x;a,b) =
    !! \int_{0}^{x} t^{a-1} (1 - t)^{b-1} dt \f$.
    interface incomplete_beta
        module procedure :: incomplete_beta_real64
        module procedure :: incomplete_beta_real32
    end interface

    !> Computes the digamma function \f$ \psi(x) = 
    !! \frac{d}{dx}\left( \ln \left( \Gamma \left( x \right) \right) 
    !! \right) \f$.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! real(real64) function digamma(real(real64) x)
    !! real(real32) function digamma(real(real32) x)
    !! @endcode
    !!
    !! @param[in] x The value at which to evaluate the function.
    !! @return The function value.
    interface digamma
        module procedure :: digamma_real64
        module procedure :: digamma_real32
    end interface
    
    !> Computes the "upper" incomplete gamma function 
    !! \f$ \Gamma(a, x) = \int_{x}^{\infty} t^{a-1} e^{-t} \,dt \f$.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! real(real64) incomplete_gamma_upper(real(real64) a, real(real64) x)
    !! real(real32) incomplete_gamma_upper(real(real32) a, real(real32) x)
    !! @endcode
    !!
    !! @param[in] a The coefficient value.
    !! @param[in] x The value at which to evaluate the function.
    !! @return The function value.
    interface incomplete_gamma_upper
        module procedure :: incomplete_gamma_upper_real64
        module procedure :: incomplete_gamma_upper_real32
    end interface

    !> Computes the "lower" incomplete gamma function 
    !! \f$ \gamma(a, x) = \int_{0}^{x} t^{a-1} e^{-t} \,dt \f$.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! real(real64) incomplete_gamma_lower(real(real64) a, real(real64) x)
    !! real(real32) incomplete_gamma_lower(real(real32) a, real(real32) x)
    !! @endcode
    !!
    !! @param[in] a The coefficient value.
    !! @param[in] x The value at which to evaluate the function.
    !! @return The function value.
    interface incomplete_gamma_lower
        module procedure :: incomplete_gamma_lower_real64
        module procedure :: incomplete_gamma_lower_real32
    end interface

    ! special_functions_beta.f90
    interface
        pure elemental module function beta_real64(a, b) result(rst)
            real(real64), intent(in) :: a, b
            real(real64) :: rst
        end function

        pure elemental module function beta_real32(a, b) result(rst)
            real(real32), intent(in) :: a, b
            real(real32) :: rst
        end function

        pure elemental module function regularized_beta_real64(a, b, x) result(rst)
            real(real64), intent(in) :: a, b, x
            real(real64) :: rst
        end function

        pure elemental module function regularized_beta_real32(a, b, x) result(rst)
            real(real32), intent(in) :: a, b, x
            real(real32) :: rst
        end function

        pure elemental module function incomplete_beta_real64(a, b, x) result(rst)
            real(real64), intent(in) :: a, b, x
            real(real64) :: rst
        end function

        pure elemental module function incomplete_beta_real32(a, b, x) result(rst)
            real(real32), intent(in) :: a, b, x
            real(real32) :: rst
        end function
    end interface

    ! special_functions_digamma.f90
    interface
        pure elemental module function digamma_real64(x) result(rst)
            real(real64), intent(in) :: x
            real(real64) :: rst
        end function

        pure elemental module function digamma_real32(x) result(rst)
            real(real32), intent(in) :: x
            real(real32) :: rst
        end function
    end interface

    ! special_functions_gamma.f90
    interface
        pure elemental module function incomplete_gamma_upper_real64(a, x) &
            result(rst)
            real(real64), intent(in) :: a, x
            real(real64) :: rst
        end function

        pure elemental module function incomplete_gamma_upper_real32(a, x) &
            result(rst)
            real(real32), intent(in) :: a, x
            real(real32) :: rst
        end function

        pure elemental module function incomplete_gamma_lower_real64(a, x) &
            result(rst)
            real(real64), intent(in) :: a, x
            real(real64) :: rst
        end function

        pure elemental module function incomplete_gamma_lower_real32(a, x) &
            result(rst)
            real(real32), intent(in) :: a, x
            real(real32) :: rst
        end function
    end interface

! ------------------------------------------------------------------------------
end module
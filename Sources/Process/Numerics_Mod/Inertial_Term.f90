!==============================================================================!
  subroutine Numerics_Mod_Inertial_Term(phi, coef, A, b, dt)
!------------------------------------------------------------------------------!
!>  The subroutine Numerics_Mod_Inertial_Term within the Numerics_Mod module
!>  discretizes the inertial term in conservation equations. It computes the
!>  contributions of the inertial term to the coefficients of the linear system
!>  being solved and adjusts the source term accordingly.
!------------------------------------------------------------------------------!
!   Functionality:                                                             !
!                                                                              !
!   * Grid aliasing: Establishes a connection to the grid associated with the  !
!     variable 'phi' for spatial discretization and volume calculations.       !
!   * Time discretization: Incorporates time discretization schemes such as    !
!     LINEAR and PARABOLIC, adapting the discretization to the chosen          !
!     temporal integration method.                                             !
!   * Coefficient adjustment: Modifies the diagonal coefficients of the        !
!     system matrix 'A' to account for the inertial contributions, enhancing   !
!     the accuracy of time-dependent simulations.                              !
!   * Source term modification: Updates the source term 'b' based on the       !
!     inertial effects, ensuring the consistency of the discretized equations. !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Arguments]--------------------------------!
  type(Var_Type)    :: phi  !! transported variable
  real              :: coef(-phi % pnt_grid % n_bnd_cells:  &
                             phi % pnt_grid % n_cells)
                            !! transport coefficient, like density for momentum,
                            !! or density times thermal capacity for enthalpy
  type(Matrix_Type) :: A    !! system matrix
  real              :: b(:) !! right-hand side, source term
  real, intent(in)  :: dt   !! time step
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  real                     :: a0
  integer                  :: c
!==============================================================================!

  ! Take alias to Grid
  Grid => phi % pnt_grid

  ! Two time levels; Linear interpolation
  if(phi % td_scheme .eq. LINEAR) then
    do c = Cells_In_Domain()
      a0 = coef(c) * Grid % vol(c) / dt
      A % val(A % dia(c)) = A % val(A % dia(c)) + a0
      b(c)  = b(c) + a0 * phi % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(phi % td_scheme .eq. PARABOLIC) then
    do c = Cells_In_Domain()
      a0 = coef(c) * Grid % vol(c) / dt
      A % val(A % dia(c)) = A % val(A % dia(c)) + 1.5 * a0
      b(c)  = b(c) + 2.0 * a0 * phi % o(c) - 0.5 * a0 * phi % oo(c)
    end do
  end if

  end subroutine

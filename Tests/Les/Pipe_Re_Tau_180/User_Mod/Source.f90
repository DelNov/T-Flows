!==============================================================================!
  subroutine User_Mod_Source(flow, phi, a_matrix, b_vector)
!------------------------------------------------------------------------------!
!   This is a prototype of a function for customized source for scalar.        !
!   It is called from "Compute_Scalar" function, just before calling the       !
!   linear solver.  Both system matrix ("a_matrix") and right hand side        !
!   vector ("b_vector") are sent should the user want to stabilize the         !
!   system for always positive variables, for example.                         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod,   only: Grid_Type
  use Field_Mod
  use Comm_Mod
  use Var_Mod,    only: Var_Type
  use Bulk_Mod,   only: Bulk_Type
  use Const_Mod,  only: PI
  use Matrix_Mod, only: Matrix_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  type(Var_Type)           :: phi
  type(Matrix_Type)        :: a_matrix
  real, dimension(:)       :: b_vector
!----------------------------------[Locals]------------------------------------!
  type(Var_Type),  pointer :: u, v, w, t
  type(Grid_Type), pointer :: grid
  type(Bulk_Type), pointer :: bulk
  integer                  :: c
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  u    => flow % u
  v    => flow % v
  w    => flow % w
  t    => flow % t
  bulk => flow % bulk

  !-----------------------------------------------------! 
  !                                                     !
  !   Set source depending on which variable you have   !
  !                                                     !
  !-----------------------------------------------------! 
  !
  !
  !   /
  !  |
  !  | dT/dz * Uz * volume
  !  |
  ! /
  !  dT/dz is derived from condition that there is no
  !  change of energy in the system. It means
  !  mass_flux * cp * dT - Q * Area = 0, Area = 2 r pi * dz ,
  !   dT/dz = 2 r pi Q / (mass_flux*cp). 
  !  where Q is heat flux through the wall.

  !-------------------------------!
  !  Set source for temperature   !
  !-------------------------------!
  if( phi % name .eq. 'T' ) then  

    call Comm_Mod_Global_Sum_Real(heat_flux)
    call Comm_Mod_Global_Sum_Real(heated_area)
    heat_flux = heat_flux / (heated_area + TINY)
    heat      = heat_flux * heated_area
     
    do c = 1, grid % n_cells
      b_vector(c) = b_vector(c)   &
                  - 2.0 * pi * heat_flux * w % n(c) &
      / (bulk % flux_z * capacity) * grid % vol(c)

    end do
  end if

  end subroutine  ! fourth level comments

!==============================================================================!
  subroutine Var_Mod_Allocate_Vol_Flux(name_phi, name_flux, phi, grid)
!------------------------------------------------------------------------------!
!   This is to allocate a variable for a solution with usual algorithm.        !
!   Variables such as velocities and pressures should be allocated with it.    !
!                                                                              !
!   One could think of storing pointer to the grid as well.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*)        :: name_phi
  character(len=*)        :: name_flux
  type(Var_Type)          :: phi
  type(Grid_Type), target :: grid
!==============================================================================!

  ! Store variable name
  phi % name      = name_phi
  phi % flux_name = name_flux

  ! Store grid for which the variable is defined
  phi % pnt_grid => grid

  ! Values new (n), last iteration (o), at old time (oo)
  allocate (phi % n (grid % n_faces));  phi % n  = 0.
  allocate (phi % o (grid % n_faces));  phi % o  = 0.
  allocate (phi % oo(grid % n_faces));  phi % oo = 0.

  end subroutine

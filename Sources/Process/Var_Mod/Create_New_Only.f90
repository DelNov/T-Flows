!==============================================================================!
  subroutine Var_Mod_Create_New_Only(phi, Grid, name_phi)
!------------------------------------------------------------------------------!
!   This is to allocate a simplified uknown, holding only current value,       !
!   such as pressure for example.                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type)          :: phi
  type(Grid_Type), target :: Grid
  character(len=*)        :: name_phi
!==============================================================================!

  ! Store Grid for which the variable is defined
  phi % pnt_grid => Grid

  ! Store variable name
  phi % name = name_phi

  ! Values in the new (n) time step
  allocate (phi % n(-Grid % n_bnd_cells:Grid % n_cells));  phi % n = 0.0

  ! Gradients
  allocate (phi % x(-Grid % n_bnd_cells:Grid % n_cells));  phi % x = 0.0
  allocate (phi % y(-Grid % n_bnd_cells:Grid % n_cells));  phi % y = 0.0
  allocate (phi % z(-Grid % n_bnd_cells:Grid % n_cells));  phi % z = 0.0

  ! Variable's boundary value
  allocate (phi % b(-Grid % n_bnd_cells:Grid % n_cells));  phi % b = 0.

  end subroutine

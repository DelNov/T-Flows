!==============================================================================!
  subroutine Var_Mod_Allocate_New_Only(name_phi, phi, grid)
!------------------------------------------------------------------------------!
!   This is to allocate a simplified uknown, holding only current value,       !
!   such as pressure for example.                                              !
!                                                                              !
!   One could think of storing pointer to the grid as well.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*)        :: name_phi
  type(Var_Type)          :: phi
  type(Grid_Type), target :: grid
!==============================================================================!

  ! Store variable name
  phi % name = name_phi

  ! Store grid for which the variable is defined
  phi % pnt_grid => grid

  ! Values in the new (n) time step
  allocate (phi % n (-grid % n_bnd_cells : grid % n_cells));   phi % n = 0.

  end subroutine

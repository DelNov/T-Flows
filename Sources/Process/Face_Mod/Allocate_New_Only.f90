!==============================================================================!
  subroutine Face_Mod_Allocate_New_Only(phi, grid, name_phi)
!------------------------------------------------------------------------------!
!   This is to allocate a face-based uknown, holding only current value.       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Face_Type)         :: phi
  type(Grid_Type), target :: grid
  character(len=*)        :: name_phi
!==============================================================================!

  ! Store grid for which the variable is defined
  phi % pnt_grid => grid

  ! Store variable name
  phi % name = name_phi

  ! Values in the new (n) time step
  allocate (phi % n(grid % n_faces));  phi % n = 0.0

  end subroutine


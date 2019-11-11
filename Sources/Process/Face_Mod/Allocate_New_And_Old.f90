!==============================================================================!
  subroutine Face_Mod_Allocate_New_Only(name_phi, phi, grid)
!------------------------------------------------------------------------------!
!   This is to allocate a face-based uknown, holding only current value.       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*)        :: name_phi
  type(Face_Type)         :: phi
  type(Grid_Type), target :: grid
!==============================================================================!

  ! Store variable name
  phi % name = name_phi

  ! Store grid for which the variable is defined
  phi % pnt_grid => grid

  ! Values in the new (n) time step
  allocate (phi % n(grid % n_faces));  phi % n = 0.0

  end subroutine

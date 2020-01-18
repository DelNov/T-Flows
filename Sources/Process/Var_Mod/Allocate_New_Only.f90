!==============================================================================!
  subroutine Var_Mod_Allocate_New_Only(phi, grid, name_phi)
!------------------------------------------------------------------------------!
!   This is to allocate a simplified uknown, holding only current value,       !
!   such as pressure for example.                                              !
!                                                                              !
!   One could think of storing pointer to the grid as well.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type)          :: phi
  type(Grid_Type), target :: grid
  character(len=*)        :: name_phi
!==============================================================================!

  ! Store grid for which the variable is defined
  phi % pnt_grid => grid

  ! Store variable name
  phi % name = name_phi

  ! Values in the new (n) time step
  allocate (phi % n(-grid % n_bnd_cells : grid % n_cells));  phi % n = 0.0

  ! Gradients
  allocate (phi % x(-grid % n_bnd_cells : grid % n_cells));  phi % x = 0.0
  allocate (phi % y(-grid % n_bnd_cells : grid % n_cells));  phi % y = 0.0
  allocate (phi % z(-grid % n_bnd_cells : grid % n_cells));  phi % z = 0.0

  ! Variable's boundary value
  allocate (phi % b(-grid % n_bnd_cells: -1));  phi % b = 0.

  end subroutine

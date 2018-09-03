!==============================================================================!
  subroutine Allocate_Additional(grid)
!------------------------------------------------------------------------------!
!   Allocates additional memory for Divisor                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Gen_Mod
  use Div_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!==============================================================================!

  ! Variables declared in Gen_Mod.  Here new numbers for cells, faces and nodes
  allocate (new_c(-grid % n_bnd_cells-1:grid % n_cells));  new_c = 0
  allocate (new_f( grid % n_faces));                       new_f = 0
  allocate (new_n( grid % n_nodes));                       new_n = 0

  ! Variables declared in Div_Mod
  allocate (buf_send_ind(grid % n_faces));  buf_send_ind = 0
  allocate (buf_recv_ind(grid % n_faces));  buf_recv_ind = 0
  allocate (buf_pos     (grid % n_faces));  buf_pos      = 0

  end subroutine

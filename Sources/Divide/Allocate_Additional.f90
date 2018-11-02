!==============================================================================!
  subroutine Allocate_Additional(grid)
!------------------------------------------------------------------------------!
!   Allocates additional memory for Divisor                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Div_Mod,  only: buf_send_ind, buf_recv_ind, buf_pos
  use Grid_Mod, only: Grid_Type,  &
                      Grid_Mod_Allocate_New_Numbers
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!==============================================================================!

  ! New nodes, cells and face numbers
  call Grid_Mod_Allocate_New_Numbers(grid,                &
                                     grid % n_nodes,      &
                                     grid % n_bnd_cells,  &
                                     grid % n_cells,      &
                                     grid % n_faces)

  ! Variables declared in Div_Mod
  allocate (buf_send_ind(grid % n_faces));  buf_send_ind = 0
  allocate (buf_recv_ind(grid % n_faces));  buf_recv_ind = 0
  allocate (buf_pos     (grid % n_faces));  buf_pos      = 0

  end subroutine

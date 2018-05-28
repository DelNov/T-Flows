!==============================================================================!
  subroutine Allocate_Additional(grid)
!------------------------------------------------------------------------------!
!   Allocates additional memory for Divisor                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use gen_mod 
  use Div_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, n
  character(len=80) :: dum_s, name_in
!==============================================================================!

  ! Variables declared in gen_mod.h90:
  allocate (new_c(-grid % n_bnd_cells-1:grid % n_cells));  new_c = 0 
  allocate (new_f( grid % n_faces));                       new_f = 0
  allocate (new_n( grid % n_nodes));                       new_n = 0 

  ! Variables declared in div.h90:
  allocate (ix(-grid % n_bnd_cells:grid % n_cells));   ix  = 0
  allocate (iy(-grid % n_bnd_cells:grid % n_cells));   iy  = 0
  allocate (iz(-grid % n_bnd_cells:grid % n_cells));   iz  = 0
  allocate (iin(-grid % n_bnd_cells:grid % n_cells));  iin = 0
  allocate (criter(grid % n_cells));  criter=0

  allocate (proces      (-grid % n_bnd_cells:grid % n_cells));  proces       = 0
  allocate (buf_send_ind(grid % n_faces));  buf_send_ind = 0
  allocate (buf_recv_ind(grid % n_faces));  buf_recv_ind = 0
  allocate (buf_pos     (grid % n_faces));  buf_pos      = 0

  end subroutine

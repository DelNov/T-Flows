!==============================================================================!
  subroutine Allocate_Memory(grid) 
!------------------------------------------------------------------------------!
!   Allocates memory for the convertor.                                        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use gen_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!==============================================================================!

  ! Allocate memory 
  ! =--> carefull: there is no checking!
  call Grid_Mod_Allocate_Nodes(grid, grid % n_nodes) 
  call Grid_Mod_Allocate_Cells(grid, grid % n_cells,   grid % n_bnd_cells) 
  call Grid_Mod_Allocate_Faces(grid, grid % n_cells*5) 

  allocate(grid % bnd_cond % color(-grid % n_bnd_cells:-1))
  grid % bnd_cond % color=0

  grid % n_copy = grid % n_faces  ! I believe it is n_cells * 5 at this point
  allocate(grid % bnd_cond % copy_c( -grid % n_bnd_cells:-1))
  grid % bnd_cond % copy_c = 0
  allocate(grid % bnd_cond % copy_s(2,grid % n_copy))
  grid % bnd_cond % copy_s=0

  allocate(new_n( grid % n_nodes));                      new_n=0  
  allocate(new_c(-grid % n_bnd_cells-1:grid % n_cells)); new_c=0  
  allocate(new_f( grid % n_cells*5));                    new_f=0  

  end subroutine

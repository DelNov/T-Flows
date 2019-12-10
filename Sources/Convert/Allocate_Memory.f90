!==============================================================================!
  subroutine Allocate_Memory(grid, mesh_city)
!------------------------------------------------------------------------------!
!   Allocates memory for the convertor.                                        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod, only: Grid_Type,                &
                      Grid_Mod_Allocate_Nodes,  &
                      Grid_Mod_Allocate_Cells,  &
                      Grid_Mod_Allocate_Faces,  &
                      Grid_Mod_Allocate_New_Numbers
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  logical         :: mesh_city  ! are we converting a mesh for a city?
!-----------------------------------[Locals]-----------------------------------!
  integer :: f = 1
!==============================================================================!

  if(mesh_city) f = 3

  ! Allocate memory 
  ! =--> carefull: there is no checking!
  call Grid_Mod_Allocate_Nodes(grid, grid % n_nodes)
  call Grid_Mod_Allocate_Cells(grid, grid % n_cells,   grid % n_bnd_cells * f)
  call Grid_Mod_Allocate_Faces(grid, grid % n_cells*5)
  call Grid_Mod_Allocate_New_Numbers(grid,                    &
                                     grid % n_nodes,          &
                                     grid % n_bnd_cells * f,  &
                                     grid % n_cells,          &
                                     grid % n_faces)

  allocate(grid % bnd_cond % color(-grid % n_bnd_cells:-1))
  grid % bnd_cond % color = 0

  end subroutine

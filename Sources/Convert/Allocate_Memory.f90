!==============================================================================!
  subroutine Allocate_Memory(grid)
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
!-----------------------------------[Locals]-----------------------------------!
  integer, parameter :: F = 3   ! workaround to allocate more memory for bnds
!==============================================================================!

  ! Allocate memory 
  ! =--> carefull: there is no checking!
  call Grid_Mod_Allocate_Nodes(grid, grid % n_nodes)
  call Grid_Mod_Allocate_Cells(grid, grid % n_cells,   grid % n_bnd_cells * F)
  call Grid_Mod_Allocate_Faces(grid, grid % n_cells*5, 0)
  call Grid_Mod_Allocate_New_Numbers(grid,                    &
                                     grid % n_nodes,          &
                                     grid % n_bnd_cells * F,  &
                                     grid % n_cells,          &
                                     grid % n_faces)

  allocate(grid % bnd_cond % color(-grid % n_bnd_cells * F:-1))
  grid % bnd_cond % color = 0

  ! (Dirty) trick to allocate additional memory for cities
  grid % n_bnd_cells = grid % n_bnd_cells / F

  end subroutine

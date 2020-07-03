!==============================================================================!
  subroutine Allocate_Additional(grid)
!------------------------------------------------------------------------------!
!   Allocates additional memory for Divisor                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
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
                                     grid % n_faces + grid % n_shadows)

  end subroutine

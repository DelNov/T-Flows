!==============================================================================!
  subroutine Allocate_Memory(grid)
!------------------------------------------------------------------------------!
! Alocates memory for working arrays and communicator                          !
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Grid_Mod, only: Grid_Type
  use Work_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!==============================================================================!

  ! Working arrays
  call Work_Mod_Allocate_Real_Cells(grid, 30)
  call Work_Mod_Allocate_Real_Faces(grid,  1)
  call Work_Mod_Allocate_Real_Nodes(grid,  1)

  call Work_Mod_Allocate_Integer_Cells(grid, 4)

  ! Variables defined in Comm_Mod.h90:
  call Comm_Mod_Allocate_Memory(grid)

  end subroutine

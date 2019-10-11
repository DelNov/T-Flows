!==============================================================================!
  subroutine Work_Mod_Allocate(grid, rc, rf, rn, ic)
!------------------------------------------------------------------------------!
!   Alocates memory for working arrays and communicator                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer :: rc  ! number of real cell
  integer :: rf  ! number of real faces
  integer :: rn  ! number of real nodes
  integer :: ic  ! number of integer cells
!==============================================================================!

  call Work_Mod_Allocate_Real_Cells(grid, rc)
  call Work_Mod_Allocate_Real_Faces(grid, rf)
  call Work_Mod_Allocate_Real_Nodes(grid, rn)

  call Work_Mod_Allocate_Integer_Cells(grid, ic)

  end subroutine

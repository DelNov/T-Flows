!==============================================================================!
  subroutine Work_Mod_Allocate(grid, rc, rf, rn, ic, if, in)
!------------------------------------------------------------------------------!
!   Alocates memory for working arrays and communicator                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid(:)
  integer         :: rc  ! number of real cell arrays
  integer         :: rf  ! number of real face arrays
  integer         :: rn  ! number of real node arrays
  integer         :: ic  ! number of integer cell arrays
  integer         :: if  ! number of integer face arrays
  integer         :: in  ! number of integer node arrays
!==============================================================================!

  call Work_Mod_Allocate_Real_Cells(grid, rc)
  call Work_Mod_Allocate_Real_Faces(grid, rf)
  call Work_Mod_Allocate_Real_Nodes(grid, rn)

  call Work_Mod_Allocate_Integer_Cells(grid, ic)
  call Work_Mod_Allocate_Integer_Faces(grid, if)
  call Work_Mod_Allocate_Integer_Nodes(grid, in)

  end subroutine

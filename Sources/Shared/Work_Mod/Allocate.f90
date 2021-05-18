!==============================================================================!
  subroutine Work_Mod_Allocate(Grid, rc, rf, rn, ic, if, in)
!------------------------------------------------------------------------------!
!   Alocates memory for working arrays and communicator                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: Grid(:)
  integer         :: rc  ! number of real cell arrays
  integer         :: rf  ! number of real face arrays
  integer         :: rn  ! number of real node arrays
  integer         :: ic  ! number of integer cell arrays
  integer         :: if  ! number of integer face arrays
  integer         :: in  ! number of integer node arrays
!==============================================================================!

  call Work_Mod_Allocate_Real_Cells(Grid, rc)
  call Work_Mod_Allocate_Real_Faces(Grid, rf)
  call Work_Mod_Allocate_Real_Nodes(Grid, rn)

  call Work_Mod_Allocate_Integer_Cells(Grid, ic)
  call Work_Mod_Allocate_Integer_Faces(Grid, if)
  call Work_Mod_Allocate_Integer_Nodes(Grid, in)

  end subroutine

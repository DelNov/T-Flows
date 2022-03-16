!==============================================================================!
  subroutine Allocate_Work(Work, Grid, rc, rf, rn, ic, if, in)
!------------------------------------------------------------------------------!
!   Alocates memory for working arrays                                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type) :: Work
  type(Grid_Type)  :: Grid(:)
  integer          :: rc  ! number of real cell arrays
  integer          :: rf  ! number of real face arrays
  integer          :: rn  ! number of real node arrays
  integer          :: ic  ! number of integer cell arrays
  integer          :: if  ! number of integer face arrays
  integer          :: in  ! number of integer node arrays
!==============================================================================!

  call Work % Allocate_Real_Cell(Grid, rc)
  call Work % Allocate_Real_Face(Grid, rf)
  call Work % Allocate_Real_Node(Grid, rn)

  call Work % Allocate_Int_Cell(Grid, ic)
  call Work % Allocate_Int_Face(Grid, if)
  call Work % Allocate_Int_Node(Grid, in)

  end subroutine

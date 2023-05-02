!==============================================================================!
  subroutine Number_Of_Swarm_Sub_Steps(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control % Read_Int_Item('NUMBER_OF_SWARM_SUB_STEPS', 60, val, verbose)

  end subroutine

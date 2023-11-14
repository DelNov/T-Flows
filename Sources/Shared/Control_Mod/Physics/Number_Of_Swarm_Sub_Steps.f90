!==============================================================================!
  subroutine Number_Of_Swarm_Sub_Steps(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the number of swarm (particle tracking) sub (time) steps.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  integer, intent(out) :: val      !! number of swarm sub (time) steps
  logical, optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Int_Item('NUMBER_OF_SWARM_SUB_STEPS', 60, val, verbose)

  end subroutine

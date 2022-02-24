!==============================================================================!
  subroutine Control_Mod_Number_Of_Swarm_Sub_Steps(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('NUMBER_OF_SWARM_SUB_STEPS', 60, &
                                  val, verbose)

  end subroutine

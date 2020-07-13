!==============================================================================!
  subroutine Control_Mod_Swarm_Save_Interval(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('SWARM_SAVE_INTERVAL', 60, &
                                  val, verbose)

  end subroutine

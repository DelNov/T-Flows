!==============================================================================!
  subroutine Control_Mod_Factor_Fictitious_Time_Vof(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('FACTOR_FICTITIOUS_TIME_VOF',  &
                                   0.1, val, verbose)

  end subroutine

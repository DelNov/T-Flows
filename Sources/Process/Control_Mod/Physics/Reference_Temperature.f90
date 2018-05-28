!==============================================================================!
  subroutine Control_Mod_Reference_Temperature(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('REFERENCE_TEMPERATURE', 0.0,  &
                                   val, verbose)

  end subroutine

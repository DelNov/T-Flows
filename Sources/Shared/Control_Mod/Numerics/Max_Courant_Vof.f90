!==============================================================================!
  subroutine Control_Mod_Max_Courant_Vof(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('MAX_COURANT_VOF',  &
                                   0.25, val, verbose)

  end subroutine

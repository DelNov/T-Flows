!==============================================================================!
  subroutine Control_Mod_Blending_Coefficient_For_User_Scalars(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('BLENDING_COEFFICIENT_FOR_USER_SCALARS',  &
                                   1.0, val, verbose)

  end subroutine

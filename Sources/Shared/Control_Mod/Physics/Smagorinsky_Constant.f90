!==============================================================================!
  subroutine Control_Mod_Smagorinsky_Constant(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('SMAGORINSKY_CONSTANT', 0.17,  &
                                   val, verbose)

  end subroutine

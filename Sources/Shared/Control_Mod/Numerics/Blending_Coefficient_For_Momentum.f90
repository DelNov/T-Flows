!==============================================================================!
  subroutine Control_Mod_Blending_Coefficient_For_Momentum(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('BLENDING_COEFFICIENT_FOR_MOMENTUM', 1.0,  &
                                 val, verbose)

  end subroutine

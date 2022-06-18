!==============================================================================!
  subroutine Control_Mod_Roughness_Coefficient(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: val(:)
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: con
!==============================================================================!

  call Control_Mod_Read_Real_Item('ROUGHNESS_COEFFICIENT', 0.0,  &
                                   con, verbose)

  ! Set the same value everywhere
  val(:) = con

  end subroutine

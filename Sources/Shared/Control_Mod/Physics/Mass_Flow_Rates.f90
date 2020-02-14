!==============================================================================!
  subroutine Control_Mod_Mass_Flow_Rates(b_x, b_y, b_z, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: b_x, b_y, b_z
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def(3)
  real :: val(3)
!==============================================================================!

  def = 0.0

  call Control_Mod_Read_Real_Array('MASS_FLOW_RATES', 3, def,  &
                                    val, verbose)

  b_x = val(1)
  b_y = val(2)
  b_z = val(3)

  end subroutine

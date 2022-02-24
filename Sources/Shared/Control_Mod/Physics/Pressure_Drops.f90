!==============================================================================!
  subroutine Control_Mod_Pressure_Drops(p_drop_x, p_drop_y, p_drop_z, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: p_drop_x, p_drop_y, p_drop_z
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def(3)
  real :: val(3)
!==============================================================================!

  data def / 0.0, 0.0, 0.0 /

  call Control_Mod_Read_Real_Array('PRESSURE_DROPS', 3, def,  &
                                    val, verbose)
  p_drop_x = val(1)
  p_drop_y = val(2)
  p_drop_z = val(3)

  end subroutine

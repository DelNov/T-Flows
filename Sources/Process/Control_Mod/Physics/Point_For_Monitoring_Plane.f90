!==============================================================================!
  subroutine Control_Mod_Point_For_Monitoring_Plane(p_x, p_y, p_z, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: p_x, p_y, p_z
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def(3)
  real :: val(3)
!==============================================================================!

  data def / 0.0, 0.0, 0.0 /

  call Control_Mod_Read_Real_Array('POINT_FOR_MONITORING_PLANE', 3, def,  &
                                    val, verbose)

  p_x = val(1)
  p_y = val(2)
  p_z = val(3)

  end subroutine

!==============================================================================!
  subroutine Control_Mod_Gravitational_Vector(grav_x, grav_y, grav_z, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: grav_x, grav_y, grav_z
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def(3)
  real :: val(3)
!==============================================================================!

  data def / 0.0, 0.0, 0.0 /

  call Control_Mod_Read_Real_Array('GRAVITATIONAL_VECTOR', 3, def,  &
                                    val, verbose)

  grav_x = val(1)
  grav_y = val(2)
  grav_z = val(3)

  end subroutine

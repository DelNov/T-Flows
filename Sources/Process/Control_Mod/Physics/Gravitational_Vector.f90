!==============================================================================!
  subroutine Control_Mod_Gravitational_Vector(v_x, v_y, v_z, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: v_x, v_y, v_z
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def(3)
  real :: val(3)
!==============================================================================!

  data def / 0.0, 0.0, 0.0 /

  call Control_Mod_Read_Real_Array('GRAVITATIONAL_VECTOR', 3, def,  &
                                    val, verbose)

  v_x = val(1)
  v_y = val(2)
  v_z = val(3)

  end subroutine

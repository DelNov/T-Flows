!==============================================================================!
  subroutine Angular_Velocity_Vector(Control, v_x, v_y, v_z, verbose)
!------------------------------------------------------------------------------!
!>  Reads the angular velocity vector from the control file.
!>  (I am not sure this function is used at all.)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control        !! parent class
  real,   intent(out) :: v_x, v_y, v_z  !! angular velocity vector component
  logical,   optional :: verbose        !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  real :: def(3)
  real :: val(3)
!==============================================================================!

  data def / 0.0, 0.0, 0.0 /

  call Control % Read_Real_Vector('ANGULAR_VELOCITY_VECTOR', 3, def,  &
                                  val, verbose)

  v_x = val(1)
  v_y = val(2)
  v_z = val(3)

  end subroutine

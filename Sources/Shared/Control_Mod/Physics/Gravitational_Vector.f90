!==============================================================================!
  subroutine Gravitational_Vector(Control, grav_x, grav_y, grav_z, verbose)
!------------------------------------------------------------------------------!
!>  Reads components of gravitational vector.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control                 !! parent class
  real,   intent(out) :: grav_x, grav_y, grav_z  !! gravity vector component
  logical,   optional :: verbose                 !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  real :: def(3)
  real :: val(3)
!==============================================================================!

  data def / 0.0, 0.0, 0.0 /

  call Control % Read_Real_Vector('GRAVITATIONAL_VECTOR', 3, def, val, verbose)

  grav_x = val(1)
  grav_y = val(2)
  grav_z = val(3)

  end subroutine

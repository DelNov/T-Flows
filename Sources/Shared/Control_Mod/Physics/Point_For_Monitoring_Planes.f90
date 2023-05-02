!==============================================================================!
  subroutine Point_For_Monitoring_Planes(Control, b_xp, b_yp, b_zp, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: b_xp, b_yp, b_zp
  logical,   optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def(3)
  real :: val(3)
!==============================================================================!

  data def / 0.0, 0.0, 0.0 /

  call Control % Read_Real_Vector('POINT_FOR_MONITORING_PLANES', 3, def,  &
                                   val, verbose)
  b_xp = val(1)
  b_yp = val(2)
  b_zp = val(3)

  end subroutine

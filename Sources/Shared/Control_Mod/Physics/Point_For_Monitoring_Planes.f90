!==============================================================================!
  subroutine Point_For_Monitoring_Planes(Control, b_xp, b_yp, b_zp, verbose)
!------------------------------------------------------------------------------!
!>  Read coordinates of the point which defines three monitoring planes
!>  (orthogonal to x, y and z direction).  These planes are used to estimate
!>  mass flow rates through the domains.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control           !! parent class
  real,   intent(out) :: b_xp, b_yp, b_zp  !! monitoring place coordinate
  logical,   optional :: verbose           !! controls output verbosity
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

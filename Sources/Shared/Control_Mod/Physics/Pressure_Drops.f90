!==============================================================================!
  subroutine Pressure_Drops(Control, p_drop_x, p_drop_y, p_drop_z, verbose)
!------------------------------------------------------------------------------!
!>  Reads pressure drops in three Cartesian coordinate directions.  These were
!>  introduced for channel flow simulations, but could be useful in any kind
!>  of simulations with periodicity in principal flow direction.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control                      !! parent class
  real,   intent(out) :: p_drop_x, p_drop_y, p_drop_z !! pressure drop component
  logical,   optional :: verbose                      !! output verbosity
!-----------------------------------[Locals]-----------------------------------!
  real :: def(3)
  real :: val(3)
!==============================================================================!

  data def / 0.0, 0.0, 0.0 /

  call Control % Read_Real_Vector('PRESSURE_DROPS', 3, def, val, verbose)

  p_drop_x = val(1)
  p_drop_y = val(2)
  p_drop_z = val(3)

  end subroutine

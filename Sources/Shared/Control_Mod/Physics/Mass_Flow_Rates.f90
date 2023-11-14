!==============================================================================!
  subroutine Mass_Flow_Rates(Control, b_x, b_y, b_z, verbose)
!------------------------------------------------------------------------------!
!>  Reads mass flow rates in three Cartesian coordinate directions.  These were
!>  introduced for channel flow simulations, but could be useful in any kind
!>  of simulations with periodicity in principal flow direction.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control        !! parent class
  real,   intent(out) :: b_x, b_y, b_z  !! mass flow rate component
  logical,   optional :: verbose        !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  real :: def(3)
  real :: val(3)
!==============================================================================!

  def = 0.0

  call Control % Read_Real_Vector('MASS_FLOW_RATES', 3, def, val, verbose)

  b_x = val(1)
  b_y = val(2)
  b_z = val(3)

  end subroutine

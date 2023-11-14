!==============================================================================!
  subroutine Reference_Temperature(Control, t_ref, verbose)
!------------------------------------------------------------------------------!
!>  Reads the value of reference temperature from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: t_ref    !! value of reference temperature
  logical,   optional :: verbose  !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  real :: def
!==============================================================================!

  def = 0.0

  call Control % Read_Real_Item('REFERENCE_TEMPERATURE', def, t_ref, verbose)

  end subroutine

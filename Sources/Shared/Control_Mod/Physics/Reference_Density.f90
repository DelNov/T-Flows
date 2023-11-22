!==============================================================================!
  subroutine Reference_Density(Control, d_ref, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: d_ref
  logical,   optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def
!==============================================================================!

  def = 0.0

  call Control % Read_Real_Item('REFERENCE_DENSITY', def, d_ref, verbose)

  end subroutine

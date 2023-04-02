!==============================================================================!
  subroutine Scalars_Diffusivity(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('SCALARS_DIFFUSIVITY', 1.0e-6, val, verbose)

  end subroutine

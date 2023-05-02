!==============================================================================!
  subroutine Max_Courant_Vof(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('MAX_COURANT_VOF', 0.25, val, verbose)

  end subroutine

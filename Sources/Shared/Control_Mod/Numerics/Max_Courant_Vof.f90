!==============================================================================!
  subroutine Max_Courant_Vof(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads maximum Courant number for VOF.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val      !! max. Courant number for VOF
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Real_Item('MAX_COURANT_VOF', 0.25, val, verbose)

  end subroutine

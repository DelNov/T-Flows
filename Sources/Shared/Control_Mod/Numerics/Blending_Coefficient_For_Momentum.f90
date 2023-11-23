!==============================================================================!
  subroutine Blending_Coefficient_For_Momentum(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads blending coefficient for momentum from control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val      !! blending coefficient value
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Real_Item('BLENDING_COEFFICIENT_FOR_MOMENTUM', 1.0,  &
                                 val, verbose)

  end subroutine

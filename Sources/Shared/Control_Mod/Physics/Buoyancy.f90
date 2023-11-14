!==============================================================================!
  subroutine Buoyancy(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads buoyancy model from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control  !! parent class
  character(SL), intent(out) :: val      !! buoyancy (thermal or density)
  logical,       optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Char_Item('BUOYANCY', 'none', val, verbose)

  call String % To_Upper_Case(val)

  end subroutine

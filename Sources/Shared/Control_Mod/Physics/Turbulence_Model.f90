!==============================================================================!
  subroutine Turbulence_Model(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads turbulence model from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control  !! parent class
  character(SL), intent(out) :: val      !! turbulence model
  logical,       optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Char_Item('TURBULENCE_MODEL', 'none', val, verbose)
  call String % To_Upper_Case(val)

  end subroutine

!==============================================================================!
  subroutine Turbulence_Model(Control, val, verbose)
!------------------------------------------------------------------------------!
!   Reading turbulence model from the control file.                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control
  character(SL), intent(out) :: val
  logical,       optional    :: verbose
!==============================================================================!

  call Control % Read_Char_Item('TURBULENCE_MODEL', 'none', val, verbose)
  call String % To_Upper_Case(val)

  end subroutine

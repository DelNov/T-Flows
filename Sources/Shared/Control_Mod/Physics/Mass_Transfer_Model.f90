!==============================================================================!
  subroutine Mass_Transfer_Model(Control, val, verbose)
!------------------------------------------------------------------------------!
!   Reading mass transfer model from the control file.                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control
  character(SL), intent(out) :: val
  logical,       optional    :: verbose
!==============================================================================!

  call Control % Read_Char_Item('MASS_TRANSFER_MODEL', 'none', val, verbose)
  call String % To_Upper_Case(val)

  end subroutine

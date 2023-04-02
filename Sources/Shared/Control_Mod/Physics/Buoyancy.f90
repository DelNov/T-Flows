!==============================================================================!
  subroutine Buoyancy(Control, val, verbose)
!------------------------------------------------------------------------------!
!   Reading buoyancy model from the control file.                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control
  character(SL), intent(out) :: val
  logical,       optional    :: verbose
!==============================================================================!

  call Control % Read_Char_Item('BUOYANCY', 'none', val, verbose)

  call String % To_Upper_Case(val)

  end subroutine

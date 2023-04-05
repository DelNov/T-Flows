!==============================================================================!
  subroutine Turbulent_Heat_Flux_Model(Control, val, verbose)
!------------------------------------------------------------------------------!
!   Reading turbulent heat flux model from the control file.                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control
  character(SL), intent(out) :: val
  logical,       optional    :: verbose
!==============================================================================!

  call Control % Read_Char_Item('TURBULENT_HEAT_FLUX_MODEL', 'SGDH',  &
                                 val, verbose)
  call String % To_Upper_Case(val)

  end subroutine

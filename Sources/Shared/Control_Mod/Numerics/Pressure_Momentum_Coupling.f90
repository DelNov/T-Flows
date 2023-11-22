!==============================================================================!
  subroutine Pressure_Momentum_Coupling(Control, name, verbose)
!------------------------------------------------------------------------------!
!   Reading pressure-momentum model from the control file.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control
  character(SL), intent(out) :: name  ! simple or piso
  logical, optional          :: verbose
!==============================================================================!

  call Control % Read_Char_Item('PRESSURE_MOMENTUM_COUPLING',  &
                                'simple', name, verbose=verbose)
  call String % To_Upper_Case(name)

  end subroutine

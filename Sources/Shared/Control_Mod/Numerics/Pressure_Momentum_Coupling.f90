!==============================================================================!
  subroutine Pressure_Momentum_Coupling(Control, name, verbose)
!------------------------------------------------------------------------------!
!>  Reads pressure-momentum model from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control  !! parent class
  character(SL), intent(out) :: name     !! method name (simple or piso)
  logical, optional          :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Char_Item('PRESSURE_MOMENTUM_COUPLING',  &
                                'simple', name, verbose=verbose)
  call String % To_Upper_Case(name)

  end subroutine

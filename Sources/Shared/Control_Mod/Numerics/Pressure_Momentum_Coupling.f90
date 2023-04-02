!==============================================================================!
  subroutine Control_Mod_Pressure_Momentum_Coupling(name, verbose)
!------------------------------------------------------------------------------!
!   Reading pressure-momentum model from the control file.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL), intent(out) :: name  ! simple or piso
  logical, optional          :: verbose
!==============================================================================!

  call Control % Read_Char_Item('PRESSURE_MOMENTUM_COUPLING',  &
                                'simple', name, verbose=verbose)
  call String % To_Upper_Case(name)

  end subroutine

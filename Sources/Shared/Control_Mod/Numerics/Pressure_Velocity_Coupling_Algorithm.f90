!==============================================================================!
  subroutine Control_Mod_Pressure_Velocity_Coupling_Algorithm(algorithm_name,  &
                                                              verbose)
!------------------------------------------------------------------------------!
!   Reading turbulence model from the control file.                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=80) :: algorithm_name
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('PRESSURE_MOMENTUM_COUPLING',  &
                                   'none', algorithm_name, verbose = .true.)
  call To_Upper_Case(algorithm_name)

  end subroutine

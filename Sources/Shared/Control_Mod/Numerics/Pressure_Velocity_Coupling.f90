!==============================================================================!
  subroutine Control_Mod_Pressure_Velocity_Coupling(name, verbose)
!------------------------------------------------------------------------------!
!   Reading pressure-velocity model from the control file.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=80) :: name  ! name of the pressure-velocity coupling algorithm
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('PRESSURE_VELOCITY_COUPLING',  &
                                   'simple', name, verbose = .true.)
  call To_Upper_Case(name)

  end subroutine

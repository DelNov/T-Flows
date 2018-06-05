!==============================================================================!
  subroutine Control_Mod_Turbulence_Wall_Treatment(verbose)
!------------------------------------------------------------------------------!
!   Reading turbulence wall treatment from the control file                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Turbulence_Mod, only: turbulence_wall_treatment,  &
                            LOW_RE,                     &
                            HIGH_RE,                    &
                            COMPOUND
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('TURBULENCE_WALL_TREATMENT', 'high_re',  &
                                   val, verbose)
  call To_Upper_Case(val)

  if(val .eq. 'LOW_RE') then
    turbulence_wall_treatment = LOW_RE

  else if(val .eq. 'HIGH_RE') then
    turbulence_wall_treatment = HIGH_RE

  else if(val .eq. 'COMPOUND') then
    turbulence_wall_treatment = COMPOUND

  else
    print *, '# Unknown turbulence wall model: ', trim(val)
    print *, '# Exiting!'
    stop 
  end if

  end subroutine

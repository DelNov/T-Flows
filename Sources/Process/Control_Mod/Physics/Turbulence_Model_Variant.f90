!==============================================================================!
  subroutine Control_Mod_Turbulence_Model_Variant(verbose)
!------------------------------------------------------------------------------!
!   Reading turbulence model variant from the control file                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Turbulence_Mod, only: turbulence_model_variant,  &
                            NONE,                      &
                            STABILIZED                   
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('TURBULENCE_MODEL_VARIANT', 'stabilized',  &
                                   val, verbose)
  call To_Upper_Case(val)

  if     (val .eq. 'NONE') then
    turbulence_model_variant = NONE

  else if(val .eq. 'STABILIZED') then
    turbulence_model_variant = STABILIZED

  else
    print *, '# Unknown turbulence model variant: ', trim(val)
    print *, '# Exiting!'
    stop 
  end if

  end subroutine

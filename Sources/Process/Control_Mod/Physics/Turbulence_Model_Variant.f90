!==============================================================================!
  subroutine Control_Mod_Turbulence_Model_Variant(verbose)
!------------------------------------------------------------------------------!
!   Reading turbulence model variant from the control file                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod,       only: this_proc, Comm_Mod_End
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
    if(this_proc < 2) then
      print *, '# Unknown turbulence model variant: ', trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End
  end if

  end subroutine

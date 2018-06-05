!==============================================================================!
  subroutine Control_Mod_Turbulence_Model_Variant(verbose)
!------------------------------------------------------------------------------!
!   Reading turbulence model variant from the control file                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Turbulence_Mod, only: turbulence_model_variant,  &
                            NONE,                      &
                            HYBRID,                    &
                            PURE,                      &
                            URANS,                     &
                            LOW_RE,                    &
                            HIGH_RE,                   &
                            WALE,                      &
                            DYNAMIC,                   &
                            SMAGORINSKY  
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('TURBULENCE_MODEL_VARIANT', 'high_re',  &
                                   val, verbose)
  call To_Upper_Case(val)

  if     (val .eq. 'NONE') then
    turbulence_model_variant = NONE
  else if(val .eq. 'HYBRID') then
    turbulence_model_variant = HYBRID
  else if(val .eq. 'PURE') then
    turbulence_model_variant = PURE
  else if(val .eq. 'URANS') then
    turbulence_model_variant = URANS
  else if(val .eq. 'LOW_RE') then
    turbulence_model_variant = LOW_RE
  else if(val .eq. 'HIGH_RE') then
    turbulence_model_variant = HIGH_RE
  else if(val .eq. 'WALE') then
    turbulence_model_variant = WALE
  else if(val .eq. 'DYNAMIC') then
    turbulence_model_variant = DYNAMIC
  else if(val .eq. 'SMAGORINSKY') then
    turbulence_model_variant = SMAGORINSKY
  else
    print *, '# Unknown turbulence model variant: ', trim(val)
    print *, '# Exiting!'
    stop 
  end if

  end subroutine

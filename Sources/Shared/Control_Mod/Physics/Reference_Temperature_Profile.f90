!==============================================================================!
  subroutine Reference_Temperature_Profile(Control, t_ref_profile_on, verbose)
!------------------------------------------------------------------------------!
!>  Reads, from the control file, whether the Boussinesq reference
!>  temperature should be a per-cell profile captured from the field at
!>  the first time step, instead of the single REFERENCE_TEMPERATURE
!>  constant.  Meant for cases with a non-trivial background vertical
!>  stratification (e.g. penetrative convection), where a single global
!>  reference value doesn't represent the ambient state.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control            !! parent class
  logical, intent(out) :: t_ref_profile_on  !! true if a captured profile
                                             !! should be used instead of
                                             !! the constant reference value
  logical, optional    :: verbose            !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control % Read_Char_Item('REFERENCE_TEMPERATURE_PROFILE',  &
                                'no', val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    t_ref_profile_on = .true.

  else if( val .eq. 'NO' ) then
    t_ref_profile_on = .false.

  else
    call Message % Error(72,                                                &
             'Unknown state for reference temperature profile: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',                       &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine

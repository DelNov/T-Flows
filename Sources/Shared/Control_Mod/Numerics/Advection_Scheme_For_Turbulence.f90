!==============================================================================!
  subroutine Control_Mod_Advection_Scheme_For_Turbulence(scheme_name, verbose)
!------------------------------------------------------------------------------!
!   Reading turbulence model from the control file.                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL), intent(out) :: scheme_name
  logical, optional          :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('ADVECTION_SCHEME_FOR_TURBULENCE', 'smart', &
                                   scheme_name, verbose)
  call String % To_Upper_Case(scheme_name)

  end subroutine

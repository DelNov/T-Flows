!==============================================================================!
  subroutine Control_Mod_Advection_Scheme_For_Vof(scheme_name, verbose)
!------------------------------------------------------------------------------!
!   Reading advection shceme for vof from control file.                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL), intent(out) :: scheme_name
  logical, optional          :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('ADVECTION_SCHEME_FOR_VOF',  &
                                  'upwind', scheme_name, verbose)
  call To_Upper_Case(scheme_name)

  end subroutine

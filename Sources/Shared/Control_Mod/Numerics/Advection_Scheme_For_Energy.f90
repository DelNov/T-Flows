!==============================================================================!
  subroutine Control_Mod_Advection_Scheme_For_Energy(scheme_name, verbose)
!------------------------------------------------------------------------------!
!   Reading turbulence model from the control file.                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL)     :: scheme_name
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('ADVECTION_SCHEME_FOR_ENERGY', 'upwind',  &
                                   scheme_name, verbose)
  call To_Upper_Case(scheme_name)

  end subroutine

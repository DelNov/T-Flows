!==============================================================================!
  subroutine Advection_Scheme_For_Vof(Control, scheme_name, verbose)
!------------------------------------------------------------------------------!
!   Reading advection shceme for vof from control file.                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control
  character(SL), intent(out) :: scheme_name
  logical, optional          :: verbose
!==============================================================================!

  call Control % Read_Char_Item('ADVECTION_SCHEME_FOR_VOF',  &
                                'upwind', scheme_name, verbose)
  call String % To_Upper_Case(scheme_name)

  end subroutine

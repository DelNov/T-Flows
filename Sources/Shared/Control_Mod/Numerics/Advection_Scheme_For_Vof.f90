!==============================================================================!
  subroutine Advection_Scheme_For_Vof(Control, scheme_name, verbose)
!------------------------------------------------------------------------------!
!>  Reads advection scheme for VOF from control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control      !! parent class
  character(SL), intent(out) :: scheme_name  !! scheme name
  logical, optional          :: verbose      !! controls output verbosity
!==============================================================================!

  call Control % Read_Char_Item('ADVECTION_SCHEME_FOR_VOF',  &
                                'upwind', scheme_name, verbose)
  call String % To_Upper_Case(scheme_name)

  end subroutine

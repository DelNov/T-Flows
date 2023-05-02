!==============================================================================!
  subroutine Time_Integration_Scheme(Control, scheme_name, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control
  character(SL), intent(out) :: scheme_name
  logical, optional          :: verbose
!==============================================================================!

  call Control % Read_Char_Item('TIME_INTEGRATION_SCHEME', 'linear',  &
                                 scheme_name, verbose)
  call String % To_Upper_Case(scheme_name)

  end subroutine

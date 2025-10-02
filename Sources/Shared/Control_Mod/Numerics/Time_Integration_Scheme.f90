!==============================================================================!
  subroutine Time_Integration_Scheme(Control, scheme_name, verbose)
!------------------------------------------------------------------------------!
!>  Reads the time-integration scheme from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control      !! parent class
  character(SL), intent(out) :: scheme_name  !! time integration scheme name
  logical, optional          :: verbose      !! controls output verbosity
!==============================================================================!

  call Control % Read_Char_Item('TIME_INTEGRATION_SCHEME', 'linear',  &
                                 scheme_name, verbose)
  call String % To_Upper_Case(scheme_name)

  end subroutine

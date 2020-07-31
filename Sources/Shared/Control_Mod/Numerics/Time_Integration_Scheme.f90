!==============================================================================!
  subroutine Control_Mod_Time_Integration_Scheme(scheme_name, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL)     :: scheme_name
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('TIME_INTEGRATION_SCHEME', 'linear',  &
                                   scheme_name, verbose)
  call To_Upper_Case(scheme_name)

  end subroutine

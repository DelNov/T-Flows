!==============================================================================!
  subroutine Control_Mod_Distance_Function_Time_Integration_Scheme(  &
                                                         scheme_name, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=80) :: scheme_name
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('DISTANCE_FUNCTION_TIME_INTEGRATION_SCHEME', &
                                  'linear', scheme_name, verbose)
  call To_Upper_Case(scheme_name)

  end subroutine

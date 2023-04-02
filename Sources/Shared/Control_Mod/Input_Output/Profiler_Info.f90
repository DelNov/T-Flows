!==============================================================================!
  subroutine Control_Mod_Profiler_Info(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL),     intent(out) :: val
  logical, optional, intent(in)  :: verbose
!==============================================================================!

  call Control % Read_Char_Item('PROFILER_INFO', 'percents', val, verbose)
  call String % To_Upper_Case(val)

  end subroutine

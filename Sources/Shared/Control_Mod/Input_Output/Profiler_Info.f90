!==============================================================================!
  subroutine Profiler_Info(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads, from the control file, how you want profiler info, in percentage
!>  points or seconds?
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)              :: Control  !! parent class
  character(SL),       intent(out) :: val      !! percents or seconds
  logical,   optional, intent(in)  :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Char_Item('PROFILER_INFO', 'percents', val, verbose)
  call String % To_Upper_Case(val)

  end subroutine

!==============================================================================!
  subroutine Report_Volume_Balance_Stop(Flow)
!------------------------------------------------------------------------------!
!   Closes file for volume balance reporting.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type) :: Flow
!==============================================================================!

  if(First_Proc()) then
    if(Flow % report_vol_balance) close(Flow % fuvbr)
  end if

  end subroutine

!==============================================================================!
  subroutine Update_By_Rank(Cpu_Timer, i_fun)
!------------------------------------------------------------------------------!
!   Updates function's timer by her rank (number)                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Cpu_Timer_Type), target :: Cpu_Timer
  integer, intent(in)           :: i_fun
!==============================================================================!

  ! Store the last time which was recorded
  Cpu_Timer % time_prev = Cpu_Timer % time_curr

  ! Refresh the value of time_curr
  call cpu_time(Cpu_Timer % time_curr)

  if(i_fun > 0) then
    Cpu_Timer % funct_time(i_fun) = Cpu_Timer % funct_time(i_fun)  &
                                  + Cpu_Timer % time_curr          &
                                  - Cpu_Timer % time_prev
  end if

  end subroutine

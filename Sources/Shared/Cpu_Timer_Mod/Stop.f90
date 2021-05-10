!==============================================================================!
  subroutine Stop(Cpu_Timer, f_name)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Cpu_Timer_Type), target :: Cpu_Timer
  character(len=*)              :: f_name
!-----------------------------------[Locals]-----------------------------------!
  integer :: f
!==============================================================================!

  !-------------------------------------!
  !   Find which function is stopping   !
  !-------------------------------------!

  ! Browse through stored functions
  do f = 1, Cpu_Timer % n_funct
    if(f_name .eq. Cpu_Timer % funct_name(f)) then
      goto 2
    end if
  end do

  ! If here, Cpu_Timer_Start wasn't invoked for this function
  print *, 'CRITICAL ERROR in ''Cpu_Timer_End'':'
  print *, 'For function ''', trim(f_name), ''', ''Cpu_Ti' // &
           'mer_Start'' wasn''t invoked.  Exiting!'
  stop

  ! Function has been found, continue
2 continue

  !-------------------------------------------------------------!
  !   Update the time for the function which is being stopped   !
  !-------------------------------------------------------------!

  ! Store the last time which was recorded
  Cpu_Timer % time_prev = Cpu_Timer % time_curr

  ! Refresh the value of time_curr
  call cpu_time(Cpu_Timer % time_curr)

  Cpu_Timer % funct_time(f) = Cpu_Timer % funct_time(f)  &
                            + Cpu_Timer % time_curr - Cpu_Timer % time_prev

  end subroutine

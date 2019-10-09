!==============================================================================!
  subroutine Cpu_Timer_Mod_Stop(f_name)
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  character(len=*) :: f_name
  integer          :: f
!==============================================================================!

  !-------------------------------------!
  !   Find which function is stopping   !
  !-------------------------------------!

  ! Browse through stored functions
  do f = 1, n_funct
    if(f_name .eq. funct_name(f)) then
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
  time_prev = time_curr     ! store the last time which was recorded
  call cpu_time(time_curr)  ! refresh the value of time_curr
  funct_time(f) = funct_time(f) + time_curr - time_prev

  end subroutine

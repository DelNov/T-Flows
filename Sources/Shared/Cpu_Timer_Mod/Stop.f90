!==============================================================================!
  subroutine Stop(Cpu_Timer, f_name)
!------------------------------------------------------------------------------!
!   Stops a function by her name                                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Cpu_Timer_Type), target :: Cpu_Timer
  character(len=*)              :: f_name
!-----------------------------------[Locals]-----------------------------------!
  integer :: i_fun
!==============================================================================!

  !----------------------------------------------------------!
  !   Find the rank (number) of function which is stopping   !
  !----------------------------------------------------------!

  ! Browse through stored functions
  do i_fun = 1, Cpu_Timer % n_functions
    if(f_name .eq. Cpu_Timer % funct_name(i_fun)) then
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
  call Cpu_Timer % Update_By_Rank(i_fun)

  !-------------------------------------------------------!
  !   Restart the function which was previously running   !
  !-------------------------------------------------------!
  Cpu_Timer % currently_running = Cpu_Timer % previously_running

  end subroutine

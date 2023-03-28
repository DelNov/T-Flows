!==============================================================================!
  subroutine Stop(Profiler, f_name)
!------------------------------------------------------------------------------!
!   Stops a function by her name                                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Profiler_Type), target :: Profiler
  character(len=*)             :: f_name
!-----------------------------------[Locals]-----------------------------------!
  integer :: i_fun
!==============================================================================!

  !----------------------------------------------------------!
  !                                                          !
  !   Find the rank (number) of function which is stopping   !
  !                                                          !
  !----------------------------------------------------------!

  !-------------------------------------!
  !   Browse through stored functions   !
  !-------------------------------------!
  do i_fun = 1, Profiler % n_functions
    if(f_name .eq. Profiler % funct_name(i_fun)) goto 2
  end do

  !----------------------------------------------!
  !   If you happen to be here, Profiler_Start   !
  !   wasn't invoked for this function           !
  !----------------------------------------------!
  print *, 'CRITICAL ERROR in ''Profiler_End'':'
  print *, 'For function ''', trim(f_name), ''', ''Cpu_Ti' // &
           'mer_Start'' wasn''t invoked.  Exiting!'
  stop

  !------------------------------------------------------------------------!
  !                                                                        !
  !   If here, function has been found and variable i_fun holds its rank   !
  !                                                                        !
  !------------------------------------------------------------------------!
2 continue

  !-------------------------------------------------------------!
  !   Update the time for the function which is being stopped   !
  !-------------------------------------------------------------!
  call Profiler % Update_By_Rank(i_fun)

  !-------------------------------------------------------!
  !   Restart the function which was previously running   !
  !-------------------------------------------------------!
  Profiler % currently_running = Profiler % previously_running(i_fun)

  end subroutine

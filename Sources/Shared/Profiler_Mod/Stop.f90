!==============================================================================!
  subroutine Stop(Prof, f_name)
!------------------------------------------------------------------------------!
!>  Stops profiling for a specific function (or code segment). It updates the
!>  profiling data and handles the transition back to the function that was
!>  running previously.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Profiler_Type), target :: Prof    !! parent class
  character(len=*)             :: f_name  !! function (or code segment) name
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
  do i_fun = 1, Prof % n_functs
    if(f_name .eq. Prof % funct_name(i_fun)) goto 2
  end do

  !----------------------------------------------!
  !   If you happen to be here, Profiler_Start   !
  !   wasn't invoked for this function           !
  !----------------------------------------------!
  call Message % Error(72,                                              &
                       'For function or section: "' // trim(f_name) //  &
                       '", Profiler % Start was not invoked. This ' //  &
                       ' error is critical.  Exiting!',                 &
                       file=__FILE__, line=__LINE__, one_proc=.true.)

  !------------------------------------------------------------------------!
  !                                                                        !
  !   If here, function has been found and variable i_fun holds its rank   !
  !                                                                        !
  !------------------------------------------------------------------------!
2 continue

  !-------------------------------------------------------------!
  !   Update the time for the function which is being stopped   !
  !-------------------------------------------------------------!
  call Prof % Update_By_Rank(i_fun)

  !-------------------------------------------------------!
  !   Restart the function which was previously running   !
  !-------------------------------------------------------!
  Prof % curr_running = Prof % prev_running(i_fun)

  end subroutine

!==============================================================================!
  subroutine Start(Profiler, f_name)
!------------------------------------------------------------------------------!
!   This subroutine is called whenever new function is invoked                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Profiler_Type), target :: Profiler
  character(len=*)              :: f_name
!-----------------------------------[Locals]-----------------------------------!
  integer :: i_fun  ! function counter
!==============================================================================!

  !---------------------------------------------------------------------------!
  !   If profiles is called for the first time, get system clock count rate   !
  !---------------------------------------------------------------------------!
  if(Profiler % n_functions .eq. 0) then
    call system_clock(count_rate = Profiler % sys_count_rate)
  end if

  !-------------------------------------------------------!
  !                                                       !
  !   Determine the function which is currently running   !
  !        (In other words, the caller function.)         !
  !                                                       !
  !-------------------------------------------------------!

  !----------------------------------------------!
  !   Check if this function was called before   !
  !----------------------------------------------!
  do i_fun = 1, Profiler % n_functions
    if(f_name .eq. Profiler % funct_name(i_fun)) goto 1
  end do

  !-----------------------------------------------------------------!
  !   If you happen to be here, the function wasn't called before   !
  !-----------------------------------------------------------------!

  ! It wasn't called before, add it to the suite of analyzed function
  Profiler % n_functions = Profiler % n_functions + 1
  Profiler % funct_name(Profiler % n_functions) = f_name

  ! Initialize times spent in the new function
  Profiler % funct_time(Profiler % n_functions) = 0.0

  ! Store currenlty running function to i_fun
  i_fun = Profiler % n_functions

  !--------------------------------------------------------------------------!
  !                                                                          !
  !   Meeting point: you are here no matter if function was called before,   !
  !   or if it was called for the first timei.  Also important to note,      !
  !   variable i_fun stores the currenly invoked function at this point.     !
  !                                                                          !
  !--------------------------------------------------------------------------!
1 continue

  !-----------------------------------------------------!
  !   Store the function which was previously running   !
  !-----------------------------------------------------!
  Profiler % previously_running(i_fun) = Profiler % currently_running

  !-------------------------------------------------------------------!
  !   Update the timer in the function which was previously running   !
  !-------------------------------------------------------------------!
  call Profiler % Update_By_Rank(Profiler % previously_running(i_fun))

  !---------------------------------------!
  !   Set the currentl running function   !
  !---------------------------------------!
  Profiler % currently_running = i_fun

  end subroutine

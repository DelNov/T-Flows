!==============================================================================!
  subroutine Start(Prof, f_name)
!------------------------------------------------------------------------------!
!>  Initializes profiling for a function. It records the function's name and
!>  updates internal counters and time stamps.  It is called whenever new
!>  function (or code segment) starts.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Profiler_Type), target :: Prof    !! parent class
  character(len=*)             :: f_name  !! function (or code segment) name
!-----------------------------------[Locals]-----------------------------------!
  integer :: i_fun  ! function counter
!==============================================================================!

  !---------------------------------------------------------------------------!
  !   If profiles is called for the first time, get system clock count rate   !
  !---------------------------------------------------------------------------!
  if(Prof % n_functs .eq. 0) then
    call system_clock(count_rate = Prof % sys_count_rate)
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
  do i_fun = 1, Prof % n_functs
    if(f_name .eq. Prof % funct_name(i_fun)) goto 1
  end do

  !-----------------------------------------------------------------!
  !   If you happen to be here, the function wasn't called before   !
  !-----------------------------------------------------------------!

  ! It wasn't called before, add it to the suite of analyzed function
  Prof % n_functs = Prof % n_functs + 1
  Prof % funct_name(Prof % n_functs) = f_name

  ! Initialize times spent in the new function
  Prof % funct_time(Prof % n_functs) = 0.0

  ! Store currenlty running function to i_fun
  i_fun = Prof % n_functs

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
  Prof % prev_running(i_fun) = Prof % curr_running

  !-------------------------------------------------------------------!
  !   Update the timer in the function which was previously running   !
  !-------------------------------------------------------------------!
  call Prof % Update_By_Rank(Prof % prev_running(i_fun))

  !---------------------------------------!
  !   Set the currentl running function   !
  !---------------------------------------!
  Prof % curr_running = i_fun

  end subroutine

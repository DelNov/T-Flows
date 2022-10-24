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

  !-----------------------------------------------------!
  !   Store the function which was previously running   !
  !-----------------------------------------------------!
  Profiler % previously_running = Profiler % currently_running

  !----------------------------------------------!
  !   If called for the first time, get system   !
  !   clock clock rate and initial clock count   !
  !----------------------------------------------!
  if(Profiler % previously_running .eq. 0) then
    call system_clock(count_rate = Profiler % sys_count_rate)
  end if

  !-------------------------------------------------------!
  !   Determine the function which is currently running   !
  !-------------------------------------------------------!

  ! Check if this function was called before
  do i_fun = 1, Profiler % n_functions
    if(f_name .eq. Profiler % funct_name(i_fun)) then  ! found the function
      Profiler % currently_running = i_fun
      goto 1
    end if
  end do

  ! It wasn't called before, add it to the suite of analyzed function
  Profiler % n_functions = Profiler % n_functions + 1
  Profiler % funct_name(Profiler % n_functions) = f_name

  ! Initialize times spent in the new function
  Profiler % funct_time(Profiler % n_functions) = 0.0

  ! Currently running function is the new one
  Profiler % currently_running = Profiler % n_functions

1 continue

  !-------------------------------------------------------------------!
  !   Update the timer in the function which was previously running   !
  !-------------------------------------------------------------------!
  call Profiler % Update_By_Rank(Profiler % previously_running)

  end subroutine

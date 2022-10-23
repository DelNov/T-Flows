!==============================================================================!
  subroutine Start(Cpu_Timer, f_name)
!------------------------------------------------------------------------------!
!   This subroutine is called whenever new function is invoked                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Cpu_Timer_Type), target :: Cpu_Timer
  character(len=*)              :: f_name
!-----------------------------------[Locals]-----------------------------------!
  integer :: i_fun  ! function counter
!==============================================================================!

  !-----------------------------------------------------!
  !   Store the function which was previously running   !
  !-----------------------------------------------------!
  Cpu_Timer % previously_running = Cpu_Timer % currently_running

  !-------------------------------------------------------!
  !   Determine the function which is currently running   !
  !-------------------------------------------------------!

  ! Check if this function was called before
  do i_fun = 1, Cpu_Timer % n_functions
    if(f_name .eq. Cpu_Timer % funct_name(i_fun)) then  ! found the function
      Cpu_Timer % currently_running = i_fun
      goto 1
    end if
  end do

  ! It wasn't called before, add it to the suite of analyzed function
  Cpu_Timer % n_functions = Cpu_Timer % n_functions + 1
  Cpu_Timer % funct_name(Cpu_Timer % n_functions) = f_name

  ! Initialize times spent in the new function
  Cpu_Timer % funct_time(Cpu_Timer % n_functions) = 0.0

  ! Currently running function is the new one
  Cpu_Timer % currently_running = Cpu_Timer % n_functions

1 continue

  !-------------------------------------------------------------------!
  !   Update the timer in the function which was previously running   !
  !-------------------------------------------------------------------!
  call Cpu_Timer % Update_By_Rank(Cpu_Timer % previously_running)

  end subroutine

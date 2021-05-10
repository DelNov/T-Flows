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
  integer :: f
!==============================================================================!

  !----------------------------!
  !   Store the old function   !
  !----------------------------!
  Cpu_Timer % old_funct = Cpu_Timer % new_funct

  !----------------------------!
  !   Determine new function   !
  !----------------------------!

  ! Check if this function was called before
  do f = 1, Cpu_Timer % n_funct
    if(f_name .eq. Cpu_Timer % funct_name(f)) then  ! found the function
      Cpu_Timer % new_funct = f
      goto 1
    end if
  end do

  ! It wasn't called before, add it
  Cpu_Timer % n_funct = Cpu_Timer % n_funct + 1
  Cpu_Timer % funct_name(Cpu_Timer % n_funct) = f_name

  ! Initialize times spent in the new function
  Cpu_Timer % funct_time(Cpu_Timer % n_funct) = 0.0

  ! Currently running function is the new one
  Cpu_Timer % new_funct = Cpu_Timer % n_funct

1 continue

  !---------------------------------------------------------------!
  !   Update the time for the function which was running before   !
  !---------------------------------------------------------------!

  ! Store the last time which was recorded
  Cpu_Timer % time_prev = Cpu_Timer % time_curr

  ! Refresh the value of time_curr
  call cpu_time(Cpu_Timer % time_curr)

  if(Cpu_Timer % n_funct > 1) then
    Cpu_Timer % funct_time(Cpu_Timer % old_funct) =  &
    Cpu_Timer % funct_time(Cpu_Timer % old_funct) + Cpu_Timer % time_curr  &
                                                  - Cpu_Timer % time_prev
  end if

  end subroutine

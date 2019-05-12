!==============================================================================!
  subroutine Cpu_Timer_Mod_Start(f_name)
!------------------------------------------------------------------------------!
!   This subroutine is called whenever new function is invoked                 !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  character(len=*) :: f_name
  integer          :: f
!==============================================================================!

  !----------------------------!
  !   Store the old function   !
  !----------------------------!
  old_funct = new_funct

  !----------------------------!
  !   Determine new function   !
  !----------------------------!

  ! Check if this function was called before
  do f = 1, n_funct
    if(f_name .eq. funct_name(f)) then  ! found the function
      new_funct = f
      goto 1
    end if
  end do

  ! It wasn't called before, add it
  n_funct             = n_funct + 1
  funct_name(n_funct) = f_name
  funct_time(n_funct) = 0.0      ! initialize times spent in the new function
  new_funct           = n_funct  ! currently running function is the new one

1 continue

  !---------------------------------------------------------------!
  !   Update the time for the function which was running before   !
  !---------------------------------------------------------------!
  time_prev = time_curr     ! store the last time which was recorded
  call cpu_time(time_curr)  ! refresh the value of time_curr
  if(n_funct > 1) then
    funct_time(old_funct) = funct_time(old_funct) + time_curr - time_prev
  end if

  end subroutine

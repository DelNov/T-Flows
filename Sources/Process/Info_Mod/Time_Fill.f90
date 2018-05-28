!==============================================================================!
  subroutine Info_Mod_Time_Fill(n, sim_time, wall_time)
!------------------------------------------------------------------------------!
!   Fills the info box with information to be written on the screen.           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod, only: this_proc    
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: n          ! time step          
  real    :: sim_time   ! simulation time
  real    :: wall_time  ! number of seconds of wall-clock time
!-----------------------------------[Locals]-----------------------------------!
  integer :: hours, minutes, seconds
!==============================================================================!

  if (this_proc < 2) then

    ! Write time step number
    write(time_info % lines(2)(33:43), '(a11)')    'Time step :'
    write(time_info % lines(2)(45:50),  '(i6)')    n

    ! Write simulation time
    write(time_info % lines(4)(27:43),    '(a17)') 'Simulation time :'
    write(time_info % lines(4)(45:53), '(1pe9.3)') sim_time
    write(time_info % lines(4)(55:57),     '(a3)') '[s]'     

    ! Write wall-clock time
    write(time_info % lines(5)(27:43), '(a17)') 'Wall-clock time :'
    hours   = floor(  wall_time  / 3600.0 )
    minutes = floor( (wall_time - 3600.0 * hours) / 60.0)
    seconds = floor(  wall_time - 3600.0 * hours - 60.0 * minutes )
    write(time_info % lines(5)(45:47), '(i3.3)')  hours
    write(time_info % lines(5)(48:48),   '(a1)')  ':'
    write(time_info % lines(5)(49:50), '(i2.2)')  minutes
    write(time_info % lines(5)(51:51),   '(a1)')  ':'
    write(time_info % lines(5)(52:53), '(i2.2)')  seconds
    write(time_info % lines(5)(55:61),   '(a7)') '[h:m:s]'     
 
  end if

  end subroutine

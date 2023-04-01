!==============================================================================!
  subroutine Info_Mod_Time_Fill(n, sim_time)
!------------------------------------------------------------------------------!
!   Fills the info box with information to be written on the screen.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(in) :: n         ! time step
  real,    intent(in) :: sim_time  ! simulation time
!-----------------------------------[Locals]-----------------------------------!
  integer  :: hours, minutes, seconds
  real(DP) :: wall_time                ! number of seconds of wall-clock time
!==============================================================================!

  ! Update current system clock and wall time
  call system_clock(Info % clock % cur)
  wall_time = real(Info % clock % cur - Info % clock % ini)  &
            / real(Info % clock % cnt)

  if(First_Proc()) then

    ! Write time step number
    write(Info % time % line(2)(55:65), '(a11)')    'Time step :'
    write(Info % time % line(2)(67:72),  '(i6)')    n

    ! Write simulation time
    write(Info % time % line(4)(45:61),    '(a17)')  'Simulation time :'
    write(Info % time % line(4)(63:71), '(1pe9.3)')  sim_time
    write(Info % time % line(4)(73:83),     '(a3)')  '[s]'

    ! Write wall-clock time
    write(Info % time % line(5)(45:61), '(a17)') 'Wall-clock time :'
    hours   = floor(  wall_time  / 3600.0 )
    minutes = floor( (wall_time - 3600.0 * hours) / 60.0)
    seconds = floor(  wall_time - 3600.0 * hours - 60.0 * minutes )
    write(Info % time % line(5)(63:65), '(i3.3)')  hours
    write(Info % time % line(5)(66:66),   '(a1)')  ':'
    write(Info % time % line(5)(67:69), '(i2.2)')  minutes
    write(Info % time % line(5)(69:69),   '(a1)')  ':'
    write(Info % time % line(5)(70:71), '(i2.2)')  seconds
    write(Info % time % line(5)(73:83),  '(a11)') '[hhh:mm:ss]'

  end if

  end subroutine

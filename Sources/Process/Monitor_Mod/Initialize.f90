!==============================================================================!
  subroutine Monitor_Mod_Initialize(grid, restart)
!------------------------------------------------------------------------------!
!   This is to set up monitoring points.                                       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Comm_Mod
  use Name_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  logical         :: restart
!----------------------------------[Calling]-----------------------------------!
  real :: Distance
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, m, l
  real              :: curr_dist, min_dist_all
  character(len=80) :: mon_file_name
  real, allocatable :: min_dist(:)
!==============================================================================!

  ! Allocate memory accordingly
  allocate(monitor % cell(monitor % n_points))
  allocate(min_dist      (monitor % n_points))

  !--------------------------------------------!
  !   Set the names for all monitoring files   !
  !--------------------------------------------!
  mon_file_name = problem_name
  mon_file_name(len_trim(problem_name)+ 1:  &
                len_trim(problem_name)+10) = '-monit.000'
  l = len_trim(mon_file_name)

  !-------------------------------!
  !   Find the monitoring cells   !
  !-------------------------------!
  do m = 1, monitor % n_points

    min_dist(m) = HUGE

    do c = 1, grid % n_cells
      curr_dist = Distance( monitor % x(m),  &
                            monitor % y(m),  &
                            monitor % z(m),  &
                            grid % xc(c),    &
                            grid % yc(c),    &
                            grid % zc(c))
      ! Current distance is smaller than the stored one 
      if(curr_dist < min_dist(m)) then
        monitor % cell(m) = c
        min_dist(m)     = curr_dist
      end if
    end do

    ! Check if smaller distance is found on another processor
    if(n_proc > 1) then
      min_dist_all = min_dist(m)
      call Comm_Mod_Global_Min_Real(min_dist_all)

      ! If there is, erase monitoring point at this_proc 
      if(abs(min_dist_all - min_dist(m)) >= TINY) then 
        monitor % cell(m) = 0
      end if
    end if
  end do

  !----------------------------------------------!
  !   Write first line in the monitoring files   !
  !----------------------------------------------!
  do m = 1, monitor % n_points

    if(monitor % cell(m) > 0) then

      write(mon_file_name(l-2:l),'(I3.3)') m

      if(.not. restart) then
        open(10+m, file = mon_file_name)
      else
        open(10+m, file = mon_file_name, position = 'append')
      endif

      write(10+m, '(a24, 3f16.6)')         &
            '# Monitoring point:',         &
            grid % xc( monitor % cell(m) ),  &
            grid % yc( monitor % cell(m) ),  &
            grid % zc( monitor % cell(m) )

    end if

  end do

  end subroutine

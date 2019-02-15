!==============================================================================!
  subroutine Monitor_Mod_Initialize(grid, restart)
!------------------------------------------------------------------------------!
!   This is to read the control file and set up monitoring points.             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  logical         :: restart
!----------------------------------[Calling]-----------------------------------!
  real :: Distance
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, m, n, l
  real              :: curr_dist, min_dist_all
  character(len=80) :: mon_file_name
  character(len=80) :: point_name
  real, allocatable :: min_dist(:)
  real              :: xyz(3), def(3)
!==============================================================================!

  ! Read number of monitoring points from control file
  call Control_Mod_Read_Int_Item('NUMBER_OF_MONITORING_POINTS', 0, &
                                  monitor % n_points, .true.)

  if(monitor % n_points .eq. 0) return

  ! Allocate memory for points
  allocate(monitor % x(monitor % n_points))
  allocate(monitor % y(monitor % n_points))
  allocate(monitor % z(monitor % n_points))

  ! Allocate memory accordingly
  allocate(monitor % cell(monitor % n_points))
  allocate(min_dist      (monitor % n_points))

  !----------------------------------------!
  !   Read monitoring points coordinates   !
  !----------------------------------------!
  do n = 1, monitor % n_points
    write(point_name, '(a,i3.3)') 'MONITORING_POINT_', n

    def = 0.  ! don't have a better idea what to set
    call Control_Mod_Read_Real_Array(point_name, 3, def, xyz, .true.)

    monitor % x(n) = xyz(1)
    monitor % y(n) = xyz(2)
    monitor % z(n) = xyz(3)
  end do

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

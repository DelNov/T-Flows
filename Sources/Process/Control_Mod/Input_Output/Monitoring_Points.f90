!==============================================================================!
  subroutine Control_Mod_Monitoring_Points(verbose)
!------------------------------------------------------------------------------!
!   Reads monintoring points from control file.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Monitor_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  integer           :: n
  real              :: xyz(3), def(3)
  character(len=80) :: point_name
!==============================================================================!

  call Control_Mod_Read_Int_Item('NUMBER_OF_MONITORING_POINTS', 0, &
                                  monitor % n_points, verbose)

  if(monitor % n_points .eq. 0) return

  ! Allocate memory for points
  allocate(monitor % x(monitor % n_points)) 
  allocate(monitor % y(monitor % n_points)) 
  allocate(monitor % z(monitor % n_points)) 

  do n = 1, monitor % n_points
    write(point_name, '(a,i3.3)') 'MONITORING_POINT_', n

    def = 0.  ! don't have a better idea what to set
    call Control_Mod_Read_Real_Array(point_name, 3, def, xyz, verbose)

    monitor % x(n) = xyz(1)
    monitor % y(n) = xyz(2)
    monitor % z(n) = xyz(3)
  end do


  end subroutine

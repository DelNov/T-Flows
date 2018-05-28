!==============================================================================!
  subroutine Control_Mod_Monitoring_Points(n_points, x, y, z, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: n_points
  real, allocatable :: x(:), y(:), z(:)
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  integer           :: n
  real              :: xyz(3), def(3)
  character(len=80) :: point_name
!==============================================================================!

  call Control_Mod_Read_Int_Item('NUMBER_OF_MONITORING_POINTS', 0, &
                                  n_points, verbose)

  if(n_points .eq. 0) return

  ! Allocate memory for points
  allocate(x(n_points)) 
  allocate(y(n_points)) 
  allocate(z(n_points)) 
  
  do n = 1, n_points
    write(point_name, '(a,i3.3)') 'MONITORING_POINT_', n

    def = 0.  ! don't have a better idea what to set
    call Control_Mod_Read_Real_Array(point_name, 3, def, xyz, verbose)

    x(n) = xyz(1)
    y(n) = xyz(2)
    z(n) = xyz(3)
  end do


  end subroutine

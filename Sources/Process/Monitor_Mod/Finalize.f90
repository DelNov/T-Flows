!==============================================================================!
  subroutine Finalize(Monitor)
!------------------------------------------------------------------------------!
!>  This subroutine is tasked with closing all monitoring files at the end of a
!>  simulation in T-Flows. It goes through each monitoring point defined in the
!>  Monitor object and, if the point is within the domain of the current
!>  processor, it closes the associated file stream.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Monitor_Type) :: Monitor  !! parent class of Monitory_Type
!-----------------------------------[Locals]-----------------------------------!
  integer :: m
!==============================================================================!

  !--------------------------------!
  !   Close the monitoring files   !
  !--------------------------------!
  do m = 1, Monitor % n_points
    if(Monitor % cell(m) > 0) close(Monitor % file_unit(m))
  end do

  end subroutine

!==============================================================================!
  subroutine Write_4_Vars(Monitor, n, flow)
!------------------------------------------------------------------------------!
!   This is to set up Monitoring points.                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Monitor_Type) :: Monitor
  integer, intent(in) :: n
  type(Field_Type)    :: flow
!-----------------------------------[Locals]-----------------------------------!
  integer :: m
!==============================================================================!

  do m = 1, Monitor % n_points
    if(Monitor % cell(m) > 0) then
      write(Monitor % file_unit(m),'(i9,4es15.5e3)')  n,                      &
                                            flow % u % n(Monitor % cell(m)),  &
                                            flow % v % n(Monitor % cell(m)),  &
                                            flow % w % n(Monitor % cell(m)),  &
                                            flow % p % n(Monitor % cell(m))
    end if
  end do

  end subroutine

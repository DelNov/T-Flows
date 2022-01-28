!==============================================================================!
  subroutine Write_5_Vars(Monitor, n, flow)
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
      if(Flow % n_scalars .eq. 0) then
        write(Monitor % file_unit(m),'(i9,5es15.5e3)')  n,                      &
                                              flow % u % n(monitor % cell(m)),  &
                                              flow % v % n(monitor % cell(m)),  &
                                              flow % w % n(monitor % cell(m)),  &
                                              flow % p % n(monitor % cell(m)),  &
                                              flow % t % n(monitor % cell(m))
      else if(Flow % n_scalars .eq. 1) then
        write(Monitor % file_unit(m),'(i9,6es15.5e3)')  n,                      &
                                              flow % u % n(monitor % cell(m)),  &
                                              flow % v % n(monitor % cell(m)),  &
                                              flow % w % n(monitor % cell(m)),  &
                                              flow % p % n(monitor % cell(m)),  &
                                              flow % t % n(monitor % cell(m)),  &
                                              flow % scalar(1) % n(monitor % cell(m))
      else if(Flow % n_scalars .eq. 2) then
        write(Monitor % file_unit(m),'(i9,7es15.5e3)')  n,                      &
                                              flow % u % n(monitor % cell(m)),  &
                                              flow % v % n(monitor % cell(m)),  &
                                              flow % w % n(monitor % cell(m)),  &
                                              flow % p % n(monitor % cell(m)),  &
                                              flow % t % n(monitor % cell(m)),  &
                                              flow % scalar(1) % n(monitor % cell(m)), &
                                              flow % scalar(2) % n(monitor % cell(m))

      end if
    end if
  end do

  end subroutine

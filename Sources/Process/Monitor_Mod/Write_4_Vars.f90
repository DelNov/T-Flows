!==============================================================================!
  subroutine Write_4_Vars(Monitor, n, Flow)
!------------------------------------------------------------------------------!
!   This is to set up Monitoring points.                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Monitor_Type) :: Monitor
  integer, intent(in) :: n
  type(Field_Type)    :: Flow
!-----------------------------------[Locals]-----------------------------------!
  integer :: m
!==============================================================================!

  do m = 1, Monitor % n_points
    if(Monitor % cell(m) > 0) then
      if(Flow % n_scalars .eq. 0) then
        write(Monitor % file_unit(m),'(i9,4es15.5e3)')  n,                      &
                                              Flow % u % n(monitor % cell(m)),  &
                                              Flow % v % n(monitor % cell(m)),  &
                                              Flow % w % n(monitor % cell(m)),  &
                                              Flow % p % n(monitor % cell(m))
      else if(Flow % n_scalars .eq. 1) then
        write(Monitor % file_unit(m),'(i9,5es15.5e3)')  n,                      &
                                              Flow % u % n(monitor % cell(m)),  &
                                              Flow % v % n(monitor % cell(m)),  &
                                              Flow % w % n(monitor % cell(m)),  &
                                              Flow % p % n(monitor % cell(m)),  &
                                              Flow % scalar(1) % n(monitor % cell(m))
      else if(Flow % n_scalars .eq. 2) then
        write(Monitor % file_unit(m),'(i9,6es15.5e3)')  n,                      &
                                              Flow % u % n(monitor % cell(m)),  &
                                              Flow % v % n(monitor % cell(m)),  &
                                              Flow % w % n(monitor % cell(m)),  &
                                              Flow % p % n(monitor % cell(m)),  &
                                              Flow % scalar(1) % n(monitor % cell(m)), &
                                              Flow % scalar(2) % n(monitor % cell(m))

      end if  
    end if
  end do

  end subroutine

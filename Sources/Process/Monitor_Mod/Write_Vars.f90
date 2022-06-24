!==============================================================================!
  subroutine Write_Vars(Monitor, curr_dt, Flow)
!------------------------------------------------------------------------------!
!   This is to write down variavles in monitoring points.                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Monitor_Type) :: Monitor
  integer, intent(in) :: curr_dt
  type(Field_Type)    :: Flow
!-----------------------------------[Locals]-----------------------------------!
  integer :: m, c, sc, fu
!==============================================================================!

  !------------------------------------------!
  !                                          !
  !   Browse through all monitoring points   !
  !                                          !
  !------------------------------------------!
  do m = 1, Monitor % n_points

    ! Take alias
    c = Monitor % cell(m)

    !-------------------------------!
    !   Cell is in this processor   !
    !-------------------------------!
    if(c > 0) then

      ! Take alias for file unit
      fu = Monitor % file_unit(m)

      ! Velocities and pressure are always written
      write(fu,'(i9)',       advance='no')  curr_dt
      write(fu,'(es15.5e3)', advance='no')  Flow % u % n(c)
      write(fu,'(es15.5e3)', advance='no')  Flow % v % n(c)
      write(fu,'(es15.5e3)', advance='no')  Flow % w % n(c)
      write(fu,'(es15.5e3)', advance='no')  Flow % p % n(c)

      ! If heat transfer, write temperature too
      if(Flow % heat_transfer) then
        write(fu,'(es15.5e3)', advance='no')  Flow % t % n(c)
      end if

      ! If there are scalars, write them as well
      do sc = 1, Flow % n_scalars
        write(fu,'(es15.5e3)', advance='no')  Flow % scalar(sc) % n(c)
      end do

      ! End the line
      write(fu,'(a)')  ' '

    end if  ! if cell > 0
  end do    ! m, monitorin point counter

  end subroutine

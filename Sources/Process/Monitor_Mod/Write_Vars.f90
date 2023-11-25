!==============================================================================!
  subroutine Write_Vars(Monitor, Flow, curr_dt)
!------------------------------------------------------------------------------!
!> This subroutine is responsible for writing the values of various flow
!> variables at the defined monitoring points during a T-Flows simulation. It
!> iterates over each monitoring point and, if the point falls within the
!> current processor's domain, writes the current values of velocity components
!> (u, v, w), pressure (p), temperature (t), and any scalar quantities to the
!> corresponding monitor files.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Monitor_Type) :: Monitor  !! parent class of Monitory_Type
  type(Field_Type)    :: Flow     !! flow field being monitored
  integer, intent(in) :: curr_dt  !! current time step
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

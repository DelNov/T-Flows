!==============================================================================!
  subroutine Monitor_Mod_Write_4_Vars(monitor, n, flow)
!------------------------------------------------------------------------------!
!   This is to set up monitoring points.                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Monitor_Type) :: monitor
  integer            :: n
  type(Field_Type)   :: flow
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type),  pointer :: u, v, w, p
  integer                  :: m
!==============================================================================!

  do m = 1, monitor % n_points
    if(monitor % cell(m) > 0) then
      write(monitor % file_unit(m),'(i9,4es15.5e3)')  n,                      &
                                            flow % u % n(monitor % cell(m)),  &
                                            flow % v % n(monitor % cell(m)),  &
                                            flow % w % n(monitor % cell(m)),  &
                                            flow % p % n(monitor % cell(m))
    end if
  end do

  end subroutine

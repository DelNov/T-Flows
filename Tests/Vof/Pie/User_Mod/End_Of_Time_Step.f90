!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(flow, turb, mult, swarm,  &
                                       n, n_stat_t, n_stat_p, time)
!------------------------------------------------------------------------------!
!          This function is computing benchmark for rising bubble.             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
  integer                       :: n         ! current time step
  integer                       :: n_stat_t  ! 1st t.s. statistics turbulence
  integer                       :: n_stat_p  ! 1st t.s. statistics particles
  real                          :: time      ! physical time
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: vof
  integer                  :: c, last_cell, fu
  real                     :: b_volume, surface, rise_velocity
  real                     :: circularity, c_position
  real, save               :: v_in, c_position_old
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  vof  => mult % vof

  if(n .eq. 1) c_position_old = 0.0

  !-------------------!
  !   Bubble volume   !
  !-------------------!
  b_volume = 0.0
  surface = 0.0
  c_position = 0.0
  rise_velocity = 0.0
  call Field_Mod_Grad_Variable(flow, vof)

  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    b_volume = b_volume + grid % vol(c) * vof % n(c)
    if (norm2((/vof % x(c),vof % y(c),vof % z(c)/)) > 1.0) then

      surface = surface + sqrt(vof % x(c) ** 2                    &
                             + vof % y(c) ** 2                    &
                             + vof % z(c) ** 2) * grid % vol(c)
    end if
    c_position = c_position + grid % yc(c) * vof % n(c) * grid % vol(c)
    rise_velocity = rise_velocity + flow % v % n(c) * vof % n(c) * grid % vol(c)
  end do

  call Comm_Mod_Global_Sum_Real(b_volume)
  call Comm_Mod_Global_Sum_Real(surface)
  call Comm_Mod_Global_Sum_Real(c_position)
  call Comm_Mod_Global_Sum_Real(rise_velocity)
  c_position    = c_position    / b_volume
  rise_velocity = rise_velocity / b_volume
  v_in = (c_position - c_position_old) / flow % dt

  IF(N > 1) THEN
    IF(THIS_PROC < 2) THEN
      PRINT *, 'C_POSITION        = ', C_POSITION
      PRINT *, 'C_POSITION_OLD    = ', C_POSITION_OLD
      PRINT *, 'RISE_VELOCITY (1) = ', RISE_VELOCITY
      PRINT *, 'RISE_VELOCITY (2) = ', V_IN
    END IF
  END IF
  c_position_old = c_position

  ! Write to file
  if (this_proc < 2) then
    call File_Mod_Append_File_For_Writing('bench-data.dat', fu)
    ! With circularity 2D:
!    write(fu,'(5(2x,e16.10e2))') time, b_volume,                    &
!                                2.0*PI/surface*sqrt(b_volume/PI),  &
!                                c_position,                        &
!                                rise_velocity
    ! With sphericity 3D:
    write(fu,'(6(2x,e16.10e2))') time, b_volume,                         &
                                PI**(1.0/3.0)*(6.0*b_volume)**(2.0/3.0)  &
                                /surface,                                &
                                c_position,                              &
                                rise_velocity,                           &
                                v_in
    close(fu)
  end if

  end subroutine

!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Flow, Turb, Vof, Swarm,  &
                                       n, n_stat_t, n_stat_p, time)
!------------------------------------------------------------------------------!
!   This function is computing benchmark for rising bubble.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  integer                  :: n         ! current time step
  integer                  :: n_stat_t  ! 1st t.s. statistics turbulence
  integer                  :: n_stat_p  ! 1st t.s. statistics particles
  real                     :: time      ! physical time
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: fun
  integer                  :: c, last_cell, fu
  real                     :: b_volume, surface, rise_velocity,  &
                              circularity, c_position
!==============================================================================!

  ! Take aliases
  grid => Flow % pnt_grid
  fun  => Vof % fun

  !-------------------!
  !   Bubble volume   !
  !-------------------!
  b_volume = 0.0
  surface = 0.0
  c_position = 0.0
  rise_velocity = 0.0
  call Flow % Grad_Variable(fun)

  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    b_volume = b_volume + grid % vol(c) * fun % n(c)
    if (norm2((/fun % x(c),fun % y(c),fun % z(c)/)) > 1.0) then

      surface = surface + sqrt(fun % x(c) ** 2                    &
                             + fun % y(c) ** 2                    &
                             + fun % z(c) ** 2) * grid % vol(c)
    end if
    c_position = c_position + grid % zc(c) * fun % n(c) * grid % vol(c)
    rise_velocity = rise_velocity + Flow % w % n(c) * fun % n(c) * grid % vol(c)
  end do

  call Comm_Mod_Global_Sum_Real(b_volume)
  call Comm_Mod_Global_Sum_Real(surface)
  call Comm_Mod_Global_Sum_Real(c_position)
  call Comm_Mod_Global_Sum_Real(rise_velocity)

  ! Write to file
  if (this_proc < 2) then
    call File % Append_For_Writing_Ascii('bench-data.dat', fu)
    !with circularity 2D:
!    write(fu,'(5(2X,E16.10E2))') time, b_volume,                    &
!                                2.0*PI/surface*sqrt(b_volume/PI),  &
!                                c_position/b_volume,               &
!                                rise_velocity/b_volume
    !with sphericity 3D:
    write(fu,'(5(2X,E16.10E2))') time, b_volume,                         &
                                PI**(1.0/3.0)*(6.0*b_volume)**(2.0/3.0)  &
                                /surface,                                &
                                c_position/b_volume,                     &
                                rise_velocity/b_volume
    close(fu)
  end if

  end subroutine

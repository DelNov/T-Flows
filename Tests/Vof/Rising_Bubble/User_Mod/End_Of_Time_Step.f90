!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Flow, Turb, Vof, Swarm,  &
                                       n_stat_t, n_stat_p)
!------------------------------------------------------------------------------!
!   This function is computing benchmark for rising bubble.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  integer                  :: n_stat_t  ! 1st t.s. statistics turbulence
  integer                  :: n_stat_p  ! 1st t.s. statistics particles
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: fun
  integer                  :: c, last_cell, fu
  real                     :: b_volume, surface, rise_velocity,  &
                              circularity, c_position
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  fun  => Vof % fun

  !-------------------!
  !   Bubble volume   !
  !-------------------!
  b_volume = 0.0
  surface = 0.0
  c_position = 0.0
  rise_velocity = 0.0
  call Flow % Grad_Variable(fun)

  do c = Cells_In_Domain()
    b_volume = b_volume + Grid % vol(c) * fun % n(c)
    if (norm2((/fun % x(c),fun % y(c),fun % z(c)/)) > 1.0) then

      surface = surface + sqrt(fun % x(c) ** 2                    &
                             + fun % y(c) ** 2                    &
                             + fun % z(c) ** 2) * Grid % vol(c)
    end if
    c_position = c_position + Grid % zc(c) * fun % n(c) * Grid % vol(c)
    rise_velocity = rise_velocity + Flow % w % n(c) * fun % n(c) * Grid % vol(c)
  end do

  call Global % Sum_Real(b_volume)
  call Global % Sum_Real(surface)
  call Global % Sum_Real(c_position)
  call Global % Sum_Real(rise_velocity)

  ! Write to file
  if (First_Proc()) then
    call File % Append_For_Writing_Ascii('bench-data.dat', fu)

    ! With sphericity 3D
    write(fu,'(5(2X,E16.10E2))') Time % Get_Time(), b_volume,             &
                                 PI**(1.0/3.0)*(6.0*b_volume)**(2.0/3.0)  &
                                 /surface,                                &
                                 c_position/b_volume,                     &
                                 rise_velocity/b_volume
    close(fu)
  end if

  end subroutine

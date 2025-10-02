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
  integer                  :: c, fu
  real                     :: b_volume, rise_velocity, c_position
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_Grid
  fun  => Vof % fun

  ! Integrate bubble volume, current position and rise velocity over cells
  b_volume      = 0.0
  c_position    = 0.0
  rise_velocity = 0.0

  do c = Cells_In_Domain()
    b_volume      = b_volume + Grid % vol(c) * fun % n(c)
    c_position    = c_position + Grid % zc(c) * fun % n(c) * Grid % vol(c)
    rise_velocity = rise_velocity + Flow % w % n(c) * fun % n(c) * Grid % vol(c)
  end do

  call Global % Sum_Real(b_volume)
  call Global % Sum_Real(c_position)
  call Global % Sum_Real(rise_velocity)

  ! Write to file
  if (First_Proc()) then
    call File % Append_For_Writing_Ascii('benchmark.dat', fu)

    write(fu,'(4(2x,e16.10e2))') Time % Get_Time(),      &
                                 b_volume,               &
                                 c_position/b_volume,    &
                                 rise_velocity/b_volume
    close(fu)
  end if

  end subroutine

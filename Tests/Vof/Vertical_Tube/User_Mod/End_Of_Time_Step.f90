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
  integer                  :: n_stat_t  ! 1st step for turbulence statist.
  integer                  :: n_stat_p  ! 1st step for particle statistics
!--------------------------------[Locals]--------------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: fun
  integer                  :: c, last_cell
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

  do c = Cells_In_Domain()
    b_volume = b_volume + Grid % vol(c) * fun % n(c)
    surface = surface + sqrt(fun % x(c) ** 2                                  &
                           + fun % y(c) ** 2                                  &
                           + fun % z(c) ** 2) * Grid % vol(c)
    c_position = c_position + Grid % zc(c) * fun % n(c) * Grid % vol(c)
    rise_velocity = rise_velocity + Flow % w % n(c) * fun % n(c) * Grid % vol(c)
  end do

  call Global % Sum_Real(b_volume)
  call Global % Sum_Real(surface)
  call Global % Sum_Real(c_position)
  call Global % Sum_Real(rise_velocity)

  if (First_Proc()) then
!    !with circularity:
!    open(unit = 70359, file='Bench-data.dat',position='APPEND')
!      write(70359,*) b_volume, 2.0*PI/surface*sqrt(b_volume/PI),              &
!                 c_position/b_volume, rise_velocity/b_volume
!    close(70359)

    !with sphericity:
    open(unit = 70359, file='bench-data.dat',position='APPEND')
      write(70359,*) b_volume, PI**(1.0/3.0)*(6.0*b_volume)**(2.0/3.0)/surface,         &
                 c_position/b_volume, rise_velocity/b_volume
    close(70359)
  end if

  end subroutine

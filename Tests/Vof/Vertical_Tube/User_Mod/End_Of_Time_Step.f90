!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Flow, turb, Vof, swarm,  &
                                       n, n_stat_t, n_stat_p, time)
!------------------------------------------------------------------------------!
!   This function is computing benchmark for rising bubble.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: swarm
  integer                  :: n         ! time step
  integer                  :: n_stat_t  ! 1st step for turbulence statist.
  integer                  :: n_stat_p  ! 1st step for particle statistics
  real                     :: time      ! physical time
!--------------------------------[Locals]--------------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: fun
  integer                  :: c, last_cell
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

  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    b_volume = b_volume + grid % vol(c) * fun % n(c)
    surface = surface + sqrt(fun % x(c) ** 2                                  &
                           + fun % y(c) ** 2                                  &
                           + fun % z(c) ** 2) * grid % vol(c)
    c_position = c_position + grid % zc(c) * fun % n(c) * grid % vol(c)
    rise_velocity = rise_velocity + Flow % w % n(c) * fun % n(c) * grid % vol(c)
  end do

  call Comm_Mod_Global_Sum_Real(b_volume)  
  call Comm_Mod_Global_Sum_Real(surface)  
  call Comm_Mod_Global_Sum_Real(c_position)  
  call Comm_Mod_Global_Sum_Real(rise_velocity)  

  if (this_proc < 2) then
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

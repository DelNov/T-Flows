!==============================================================================!
  subroutine Compute_Benchmark(mult, dt)
!------------------------------------------------------------------------------!
!   This function aftar VOF is solved                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Field_Type),     pointer :: flow
  type(Grid_Type),      pointer :: grid
  type(Var_Type),       pointer :: vof
  real                          :: dt    ! time step
!--------------------------------[Locals]--------------------------------------!
  character(len=1)              :: charI
  integer                       :: c, last_cell
  real                          :: b_volume, surface, rise_velocity,           &
                                   circularity, c_position
!==============================================================================!
  
  ! Take aliases
  flow     => mult % pnt_flow
  grid     => flow % pnt_grid
  vof      => mult % vof

  !-------------------!
  !   Bubble volume   ! 
  !-------------------!

  b_volume = 0.0
  surface = 0.0
  c_position = 0.0
  rise_velocity = 0.0

  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    b_volume = b_volume + grid % vol(c) * vof % n(c)
    surface = surface + sqrt(vof % x(c) ** 2                                  &
                           + vof % y(c) ** 2                                  &
                           + vof % z(c) ** 2) * grid % vol(c)
    c_position = c_position + grid % zc(c) * vof % n(c) * grid % vol(c)
    rise_velocity = rise_velocity + flow % w % n(c) * vof % n(c) * grid % vol(c)
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
    open(unit = 70359, file='Bench-data.dat',position='APPEND')
      write(70359,*) b_volume, PI**(1.0/3.0)*(6.0*b_volume)**(2.0/3.0)/surface,         &
                 c_position/b_volume, rise_velocity/b_volume
    close(70359)
  end if

  end subroutine

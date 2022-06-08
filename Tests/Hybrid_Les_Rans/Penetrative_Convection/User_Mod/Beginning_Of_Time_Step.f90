!==============================================================================!
  subroutine User_Mod_Beginning_Of_Time_Step(Flow, Turb, Vof, Swarm,  &
                                             curr_dt, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step and has two parts:   !
!     1. Keep the temperature in the uppar part of the domain constant         !
!     2. Update physical properties                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Swarm_Type),    target :: Swarm
  integer, intent(in)         :: curr_dt  ! current time step
  real,    intent(in)         :: time     ! current physical time
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: u, v, w, t, phi
  integer                  :: c, i
  real                     :: wi, wip
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  call Flow % Alias_Energy(t)

  !----------------------------------------------------------------------!
  !                                                                      !
  !   1. Keep the temperature in the uppar part of the domain constant   !
  !                                                                      !
  !----------------------------------------------------------------------!
  do c = 1, Grid % n_cells
    if(grid % wall_dist(c) > 400.0) then
      t % n(c)  = -2.0 + 0.008 * Grid % wall_dist(c)
      t % o(c)  = t % n(c)
      t % oo(c) = t % n(c)
    end if
  end do

  !-----------------------------------!
  !                                   !
  !   2. Update physical properties   !
  !                                   !
  !-----------------------------------!
  do c = -Grid % n_bnd_cells, Grid % n_cells

    ! Browse through all table entries
    do i = 1, N_ITEMS - 1

      ! Did you find the right interval
      if(Flow % t % n(c) >= air_t(i) .and.  &
         Flow % t % n(c) <  air_t(i+1)) then

        ! If so, calculate interpolation factors ...
        wi  = (air_t(i+1) - Flow % t % n(c)) / (air_t(i+1) - air_t(i))
        wip = 1.0 - wi

        ! ... and interpolate physical properties
        Flow % density(c)      = wi * air_rho   (i)  + wip * air_rho   (i+1)
        Flow % viscosity(c)    = wi * air_mu    (i)  + wip * air_mu    (i+1)
        Flow % conductivity(c) = wi * air_lambda(i)  + wip * air_lambda(i+1)
        Flow % capacity(c)     = wi * air_cp    (i)  + wip * air_cp    (i+1)
      end if

    end do
  end do

  end subroutine


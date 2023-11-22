!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Flow, Turb, Vof, Swarm,  &
                                       n_stat_t, n_stat_p)
!------------------------------------------------------------------------------!
!   Append viscous forces to pressure drops.                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  integer, intent(in)      :: n_stat_t  ! 1st step for turbulence statist.
  integer, intent(in)      :: n_stat_p  ! 1st step for particle statistics
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w ! [m/s]
  real, pointer            :: visc(:) ! [kg/(m s)]
  integer                  :: s, c1, c2
  real                     :: nx, ny, nz
  real                     :: ut, vt, wt, ut_mag ! [m/s]
  real                     :: tau_wall_x, tau_wall_y, tau_wall_z,  &
                              tau_wall ! [(kg m)/s^2]
!==============================================================================!

  ! Take the alias to the Grid
  Grid => Flow % pnt_grid
  bulk => Flow % bulk
  visc => Flow % viscosity
  call Flow % Alias_Momentum(u, v, w)

  ! Initialize stresses at the walls
  tau_wall_x = 0.0
  tau_wall_y = 0.0
  tau_wall_z = 0.0

  !---------------------------------------!
  !   Browse through all boundary faces   !
  !---------------------------------------!
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    if(c2 < 0) then

      ! If wall boundary condition
      if(Grid % Bnd_Cond_Type(c2) .eq. WALL .or.  &
         Grid % Bnd_Cond_Type(c2) .eq. WALLFL) then

        ! Wall normal
        nx = Grid % sx(s) / Grid % s(s)
        ny = Grid % sy(s) / Grid % s(s)
        nz = Grid % sz(s) / Grid % s(s)

        ! Velocity tangential to the wall (u_t = u - u_n = u (1.0-n))
        ut = u % n(c1) * (1.0 - nx)
        vt = v % n(c1) * (1.0 - ny)
        wt = w % n(c1) * (1.0 - nz)

        ut_mag = sqrt(ut**2 + vt**2 + wt**2)

        ! Shear stress magnitude at the boundary cell
        !tau_wall = visc(c1) * ut_mag  &
        !       / Grid % wall_dist(c1) * Grid % s(s)

        !tau_wall_x = tau_wall_x + tau_wall * ut / ut_mag
        !tau_wall_y = tau_wall_y + tau_wall * vt / ut_mag
        !tau_wall_z = tau_wall_z + tau_wall * wt / ut_mag

        ! Summarize shear stress components at the boundary cell
        tau_wall_x = tau_wall_x + visc(c1) * ut &
                   / Grid % wall_dist(c1) * Grid % s(s)
        tau_wall_y = tau_wall_y + visc(c1) * vt &
                   / Grid % wall_dist(c1) * Grid % s(s)
        tau_wall_z = tau_wall_z + visc(c1) * wt &
                   / Grid % wall_dist(c1) * Grid % s(s)
      end if  ! if wall
    end if
  end do

  if( abs(bulk % flux_x_o) >= TINY ) then
    call Global % Sum_Real(tau_wall_x)
    bulk % p_drop_x = bulk % p_drop_x  &
                    + tau_wall_x / Grid % tot_vol  ! [kg/m^2/s^2]
  end if
  if( abs(bulk % flux_y_o) >= TINY ) then
    call Global % Sum_Real(tau_wall_y)
    bulk % p_drop_y = bulk % p_drop_y  &
                    + tau_wall_y / Grid % tot_vol  ! [kg/m^2/s^2]
  end if
  if( abs(bulk % flux_z_o) >= TINY ) then
    call Global % Sum_Real(tau_wall_z)
    bulk % p_drop_z = bulk % p_drop_z  &
                    + tau_wall_z / Grid % tot_vol  ! [kg/m^2/s^2]
  end if


  end subroutine

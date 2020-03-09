!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(flow, turb, mult, swarm, n, time)
!------------------------------------------------------------------------------!
!   Append viscous forces to pressure drops.                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
  integer                       :: n     ! time step
  real                          :: time  ! physical time
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: s, c1, c2
  real                     :: nx, ny, nz, ut, vt, wt, ut_mag, stress
  real                     :: tau_wall_x, tau_wall_y, tau_wall_z
!==============================================================================!

  ! Take the alias to the grid
  grid => flow % pnt_grid

  ! Initialize stresses at the walls
  tau_wall_x = 0.0  ! [(kg m)/s^2]
  tau_wall_y = 0.0  ! [(kg m)/s^2]
  tau_wall_z = 0.0  ! [(kg m)/s^2]

  !---------------------------------------!
  !   Browse through all boundary faces   !
  !---------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then

      ! If wall boundary condition
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

        ! Wall normal
        nx = grid % sx(s) / grid % s(s)
        ny = grid % sy(s) / grid % s(s)
        nz = grid % sz(s) / grid % s(s)

        ! Velocity tangential to the wall (u_t = u - u_n = u (1.0-n))
        ut = flow % u % n(c1) * (1.0 - nx)
        vt = flow % v % n(c1) * (1.0 - ny)
        wt = flow % w % n(c1) * (1.0 - nz)

        ut_mag = sqrt(ut**2 + vt**2 + wt ** 2)

        ! Stress at the boundary cell [(kg m)/s^2 = N]
        stress = flow % viscosity(c1) * ut_mag  &
               / grid % wall_dist(c1) * grid % s(s)

        tau_wall_x = tau_wall_x + stress * ut / ut_mag  ! [(kg m)/s^2 = N]
        tau_wall_y = tau_wall_y + stress * vt / ut_mag  ! [(kg m)/s^2 = N]
        tau_wall_z = tau_wall_z + stress * wt / ut_mag  ! [(kg m)/s^2 = N]
      end if  ! if wall
    end if
  end do

  call Comm_Mod_Global_Sum_Real(tau_wall_x)
  call Comm_Mod_Global_Sum_Real(tau_wall_y)
  call Comm_Mod_Global_Sum_Real(tau_wall_z)

  flow % bulk % p_drop_x = flow % bulk % p_drop_x  &
                         + tau_wall_x / grid % tot_vol  ! [kg/m^2/s^2]
  flow % bulk % p_drop_y = flow % bulk % p_drop_y  &
                         + tau_wall_y / grid % tot_vol  ! [kg/m^2/s^2]
  flow % bulk % p_drop_z = flow % bulk % p_drop_z  &
                         + tau_wall_z / grid % tot_vol  ! [kg/m^2/s^2]

  end subroutine

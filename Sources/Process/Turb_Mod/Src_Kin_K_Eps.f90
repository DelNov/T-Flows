!==============================================================================!
  subroutine Turb_Mod_Src_Kin_K_Eps(turb, sol)
!------------------------------------------------------------------------------!
!   Computes the source terms in kin transport equation for k-epsilon model    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type),   target :: turb
  type(Solver_Type), target :: sol
!---------------------------------[Calling]------------------------------------!
  real :: Roughness_Coefficient
  real :: Tau_Wall_Low_Re
  real :: Tau_Wall_Rough_Walls
  real :: Y_Plus_Low_Re
  real :: Y_Plus_Rough_Walls
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: flow
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, v, w
  type(Var_Type),    pointer :: kin, eps, ut, vt, wt
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: c, c1, c2, s
  real                       :: u_tan, u_tau
  real                       :: kin_vis  ! [m^2/s]
  real                       :: p_kin_int, p_kin_wf
  real                       :: z_o
!==============================================================================!
!   Dimensions:                                                                !
!                                                                              !
!   production    p_kin    [m^2/s^3]   | rate-of-strain  shear     [1/s]       !
!   dissipation   eps % n  [m^2/s^3]   | turb. visc.     vis_t     [kg/(m*s)]  !
!   wall shear s. tau_wall [kg/(m*s^2)]| dyn visc.       viscosity [kg/(m*s)]  !
!   density       density  [kg/m^3]    | turb. kin en.   kin % n   [m^2/s^2]   !
!   cell volume   vol      [m^3]       | length          lf        [m]         !
!   left hand s.  a        [kg/s]      | right hand s.   b         [kg*m^2/s^3]!
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!------------------------------------------------------------------------------!

  ! Take aliases
  flow => turb % pnt_flow
  grid => flow % pnt_grid
  call Field_Mod_Alias_Momentum  (flow, u, v, w)
  call Turb_Mod_Alias_K_Eps      (turb, kin, eps)
  call Turb_Mod_Alias_Heat_Fluxes(turb, ut, vt, wt)
  call Solver_Mod_Alias_System   (sol, a, b)

  !-----------------------------------------!
  !   Compute the sources in the interior   !
  !-----------------------------------------!
  ! Production source:
  do c = 1, grid % n_cells
    turb % p_kin(c) = turb % vis_t(c) * flow % shear(c)**2
    b(c) = b(c) + turb % p_kin(c) * grid % vol(c)

    a % val(a % dia(c)) = a % val(a % dia(c)) + &
         flow % density(c) * eps % n(c)/(kin % n(c) + TINY) * grid % vol(c)

    if (buoyancy) then
      turb % g_buoy(c) = -flow % beta           &
                       * (grav_x * ut % n(c) +  &
                          grav_y * vt % n(c) +  &
                          grav_z * wt % n(c))   &
                       * flow % density(c)
      b(c) = b(c) + max(0.0, turb % g_buoy(c) * grid % vol(c))
      a % val(a % dia(c)) = a % val(a % dia(c))        &
                          + max(0.0,-turb % g_buoy(c)  &
                          * grid % vol(c) / (kin % n(c) + TINY))
    end if
  end do

  !-----------------------------------------------!
  !  Compute the sources in the near wall cells   !
  !-----------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then
      kin_vis = flow % viscosity(c1) / flow % density(c1)
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

        ! Compute tangential velocity component
        u_tan = Field_Mod_U_Tan(flow, s)

        u_tau = c_mu25 * sqrt(kin % n(c1))
        turb % y_plus(c1) = Y_Plus_Low_Re(turb,                  &
                                          u_tau,                 &
                                          grid % wall_dist(c1),  &
                                          kin_vis)

        turb % tau_wall(c1) = Tau_Wall_Low_Re(turb,               &
                                              flow % density(c1), &
                                              u_tau,              &
                                              u_tan,              &
                                              turb % y_plus(c1))

        p_kin_wf  = turb % tau_wall(c1) * c_mu25 * sqrt(kin % n(c1))  &
                  / (grid % wall_dist(c1) * kappa)

        p_kin_int = turb % vis_t(c1) * flow % shear(c1)**2

        turb % p_kin(c1) = p_kin_wf

        if(turb % rough_walls) then
          z_o = Roughness_Coefficient(turb, turb % z_o_f(c1))

          turb % y_plus(c1) = Y_Plus_Rough_Walls(turb,                  &
                                                 u_tau,                 &
                                                 grid % wall_dist(c1),  &
                                                 kin_vis)

          turb % tau_wall(c1) = Tau_Wall_Rough_Walls(turb,                  &
                                                     flow % density(c1),    &
                                                     u_tau,                 &
                                                     u_tan,                 &
                                                     grid % wall_dist(c1),  &
                                                     z_o)

          turb % p_kin(c1) = turb % tau_wall(c1) * c_mu25 * sqrt(kin % n(c1))  &
                           / (kappa*(grid % wall_dist(c1) + z_o))

        end if  ! rough_walls

        b(c1) = b(c1)                                                        &
              + (turb % p_kin(c1) - turb % vis_t(c1) * flow % shear(c1)**2)  &
              * grid % vol(c1)

      end if    ! Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL or WALLFL
    end if      ! c2 < 0
  end do

  call Grid_Mod_Exchange_Real(grid, kin % n)

  end subroutine

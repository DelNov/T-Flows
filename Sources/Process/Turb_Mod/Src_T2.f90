!==============================================================================!
  subroutine Turb_Mod_Src_T2(turb, sol)
!------------------------------------------------------------------------------!
!   Computes the source terms in t2 transport equation for k-eps_t2 model      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type),   target :: turb
  type(Solver_Type), target :: sol
!---------------------------------[Calling]------------------------------------!
  real :: Y_Plus_Low_Re
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: flow
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, v, w, t
  type(Var_Type),    pointer :: kin, eps, ut, vt, wt, t2
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: c, c1, c2, s
  real                       :: kin_vis, p_t2_wall, ebf, u_tau
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
  call Field_Mod_Alias_Energy    (flow, t)
  call Turb_Mod_Alias_K_Eps      (turb, kin, eps)
  call Turb_Mod_Alias_Heat_Fluxes(turb, ut, vt, wt)
  call Solver_Mod_Alias_System   (sol, a, b)
  call Turb_Mod_Alias_T2         (turb, t2)

  !-----------------------------------------!
  !   Compute the sources in all the cells  !
  !-----------------------------------------!

  ! Temperature gradients are needed
  call Field_Mod_Grad_Variable(flow, t)

  ! Production source:
  do c = 1, grid % n_cells

    turb % p_t2(c) = - 2.0 * (  ut % n(c) * t % x(c)   &
                              + vt % n(c) * t % y(c)   &
                              + wt % n(c) * t % z(c))

    b(c) = b(c) + turb % p_t2(c) * grid % vol(c)

   ! Negative contribution
   a % val(a % dia(c)) = a % val(a % dia(c)) +  &
         2.0 * flow % density(c)  * eps % n(c)  &
             / (kin % n(c) + TINY) * grid % vol(c)

  end do

  ! Implementation of wall function approach for buoyancy-driven flows

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or. &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

        ! Kinematic viscosities
        kin_vis = flow % viscosity(c1) / flow % density(c1)

        u_tau = c_mu25 * sqrt(kin % n(c1))

        turb % y_plus(c1) = Y_Plus_Low_Re(turb, u_tau,      &
                     grid % wall_dist(c1), kin_vis)

        ebf = Turb_Mod_Ebf_Momentum(turb, c1)

        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) &
        t % q(c2) = abs(turb % con_w(c1)*(t % n(c1) &
                    - t % n(c2))/grid % wall_dist(c1))

        p_t2_wall  = t % q(c2)*c_mu_theta5*sqrt(abs(t2 % n(c1))) &
                     /(kappa_theta*c_mu25*grid % wall_dist(c1))

        b(c1) = b(c1) - turb % p_t2(c1) * grid % vol(c1)

        if(turb % y_plus(c1) > 11.0) then
          turb % p_t2(c1) = p_t2_wall
        else
          turb % p_t2(c1) = (  turb % p_t2(c1) * exp(-1.0 * ebf)  &
                             + p_t2_wall * exp(-1.0/ebf))
        end if

        b(c1) = b(c1) + turb % p_t2(c1) * grid % vol(c1)

        t2 % n(c2) = 0.0

      end if  ! Grid_Mod_Bnd_Cond_Type(grid,c2).eq.WALL or WALLFL
    end if    ! c2 < 0
  end do


  end subroutine


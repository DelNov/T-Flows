!==============================================================================!
  subroutine Turb_Mod_Vis_T_Smagorinsky(turb)
!------------------------------------------------------------------------------!
!   Calculates SGS stresses and turbulent viscosity for 'LES'.                 !
!------------------------------------------------------------------------------!
!   near(c) is the number of corresponding cell on the nearest wall.           !
!   In case that, in parallel executions, the subdomain does not have          !
!   any nearwall cells, the nearest_wall_cell(c) is zero.                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target :: turb
!---------------------------------[Calling]------------------------------------!
  real :: Roughness_Coefficient
  real :: U_Plus_Rough_Walls 
  real :: Y_Plus_Rough_Walls
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w, t
  integer                   :: c, s, c1, c2
  real                      :: nx, ny, nz
  real                      :: cs, lf, u_tau_l, u_f, nc2
  real                      :: u_tan, a_pow, b_pow, nu, dely
  real                      :: beta, pr, ebf, u_plus, pr_t, sc, z_o
!==============================================================================!

  ! Take aliases
  Flow => turb % pnt_flow
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)
  t    => Flow % t

  !---------------!
  !               !
  !   SGS terms   !
  !               !
  !---------------!
  if(Flow % buoyancy .eq. THERMALLY_DRIVEN) then
    call Flow % Grad_Variable(t)
  end if

  if(turb % model .eq. LES_SMAGORINSKY) then
    do c = 1, Grid % n_cells
      lf = Grid % vol(c)**ONE_THIRD

      ! if(nearest_wall_cell(c) .ne. 0) is needed for parallel version
      ! since the subdomains which do not "touch" wall
      ! has nearest_wall_cell(c) = 0. 
      if(turb % nearest_wall_cell(c) .ne. 0) then
        u_f = sqrt( Flow % viscosity(c)                                &
                    * sqrt(  u % n(turb % nearest_wall_cell(c)) ** 2   &
                           + v % n(turb % nearest_wall_cell(c)) ** 2   &
                           + w % n(turb % nearest_wall_cell(c)) ** 2)  &
                   / (Grid % wall_dist(turb % nearest_wall_cell(c))+TINY) )
        turb % y_plus(c) = Grid % wall_dist(c) * u_f / Flow % viscosity(c)
        cs = c_smag * (1.0 - exp(-turb % y_plus(c) / 25.0))
      else  
        cs = c_smag
      end if
      turb % vis_t(c) = Flow % density(c)  &
                      * (lf*lf)            &  ! delta^2
                      * (cs*cs)            &  ! cs^2
                      * Flow % shear(c)
    end do

  else if(turb % model .eq. LES_DYNAMIC) then
    do c = 1, Grid % n_cells
      lf = Grid % vol(c) ** ONE_THIRD
      turb % vis_t(c) = Flow % density(c)  &
                      * (lf*lf)            &  ! delta^2
                      * turb % c_dyn(c)    &  ! c_dynamic
                      * Flow % shear(c)
    end do
  else if(turb % model .eq. LES_WALE) then
    do c = 1, Grid % n_cells
      lf = Grid % vol(c)**ONE_THIRD
      turb % vis_t(c) = Flow % density(c)  &
                      * (lf*lf)            &  ! delta^2
                      * (0.5*0.5)          &  ! cs^2
                      * turb % wale_v(c)
    end do
  end if

  if(Flow % buoyancy .eq. THERMALLY_DRIVEN) then
    do c = 1, Grid % n_cells
      nc2 = - Flow % beta * (  grav_x * t % x(c)   &   
                             + grav_y * t % y(c)   &   
                             + grav_z * t % z(c))  
      turb % vis_t(c) = turb % vis_t(c)            &
             * max(1 - 2.5 * nc2 / (Flow % shear(c) + TINY), 0.0)
    end do
  end if

  !-------------------!
  !   Wall function   !
  !-------------------+--------------!
  !   Law of the wall:               !
  !                                  !
  !   u+ = yz+  for z+ < 11.81       !
  !                                  !
  !   and                            !
  !                                  !
  !   u+ = A(y+)^B   for y+ > 11.81  !
  !                                  !
  !   with: A = 8.3 and B = 1/7      !
  !                                  !
  !----------------------------------+----------!
  !   The procedure below should be activated   !
  !   only if wall function approach is used.   !
  !----------------.----------------------------! 
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    if(c2  < 0) then

      nx = Grid % sx(s) / Grid % s(s)
      ny = Grid % sy(s) / Grid % s(s)
      nz = Grid % sz(s) / Grid % s(s)

      if(Grid % Bnd_Cond_Type(c2) .eq. WALL .or.  &
         Grid % Bnd_Cond_Type(c2) .eq. WALLFL) then

        u_tan = Flow % U_Tan(s)

        a_pow = 8.3
        b_pow = 1.0/7.0
        nu = Flow % viscosity(c1) / Flow % density(c1)
        dely = Grid % wall_dist(c1)

        ! Calculate u_tau_l
        u_tau_l = ( u_tan/a_pow * (nu/dely)**b_pow ) ** (1.0/(1.0+b_pow))

        ! Calculate tau_wall (but it is never used)
        turb % tau_wall(c1) = Flow % viscosity(c1) * u_tan / dely

        ! Calculate y+
        turb % y_plus(c1)  = dely * u_tau_l / nu

        if(turb % y_plus(c1)  >=  11.81) then

          if(turb % rough_walls) then
            z_o = Roughness_Coefficient(turb, turb % z_o_f(c1))
            z_o = max(Grid % wall_dist(c1)   &   
                / (e_log * max(turb % y_plus(c1),1.0)), z_o)

            turb % y_plus(c1) = Y_Plus_Rough_Walls(turb,                  &
                                                   u_tau_l,               &
                                                   Grid % wall_dist(c1),  &
                                                   nu)
            u_plus     = U_Plus_Rough_Walls(turb, Grid % wall_dist(c1))
          end if 


          ! This one is effective viscosity
          turb % vis_w(c1) = Flow % density(c1) * u_tau_l * u_tau_l * dely  &
                           / abs(u_tan)
        else
          turb % vis_w(c1) = Flow % viscosity(c1)                &
                        +      Grid % fw(s)  * turb % vis_t(c1)  &
                        + (1.0-Grid % fw(s)) * turb % vis_t(c2)
        end if

        if(Flow % heat_transfer) then
          u_plus = u_tan / u_tau_l 


          pr_t = Turb_Mod_Prandtl_Number(turb, c1)
          pr   = Flow % Prandtl_Number(c1)          ! laminar Prandtl number
          beta = 9.24 * ((pr/pr_t)**0.75 - 1.0)     &
               * (1.0 + 0.28 * exp(-0.007*pr/pr_t))

          ! According to Toparlar et al. 2019 paper
          ! "CFD simulation of the near-neutral atmospheric boundary layer: New
          ! temperature inlet profile consistent with wall functions"

          if(turb % rough_walls) then
            beta = 0.0    
            u_plus = U_Plus_Rough_Walls(turb, Grid % wall_dist(c1))
          end if

          ebf = Turb_Mod_Ebf_Scalar(turb, c1, pr) 
          turb % con_w(c1) =    turb % y_plus(c1)                         &
                              * Flow % viscosity(c1)                      &
                              * Flow % capacity(c1)                       &
                      / (  turb % y_plus(c1) * pr * exp(-1.0 * ebf)       &
                         + (u_plus + beta) * pr_t * exp(-1.0 / ebf) + TINY)
        end if

        if(Flow % n_scalars > 0) then
          u_plus = u_tan / u_tau_l 
 
          sc   = Flow % Schmidt_Number(c1)          ! laminar Schmidt number
          beta = 9.24 * ((sc/sc_t)**0.75 - 1.0)                   &
               * (1.0 + 0.28 * exp(-0.007*sc/sc_t))

          ! According to Toparlar et al. 2019 paper
          ! "CFD simulation of the near-neutral atmospheric boundary layer: New
          ! temperature inlet profile consistent with wall functions"

          if(turb % rough_walls) then
            beta = 0.0    
            u_plus = U_Plus_Rough_Walls(turb, Grid % wall_dist(c1))
          end if

          ebf = 0.01 * (sc * turb % y_plus(c1)**4                 &
              / ((1.0 + 5.0 * sc**3 * turb % y_plus(c1)) + TINY))
          turb % diff_w(c1) =  turb % y_plus(c1)                  &
              * (Flow % viscosity(c1)/Flow % density(c1))         &
              / (turb % y_plus(c1) * sc * exp(-1.0 * ebf)         &
              + (u_plus + beta) * sc_t * exp(-1.0 / ebf) + TINY)
        end if

      end if  ! Grid % Bnd_Cond_Type(c2) .eq. WALL or WALLFL
    end if    ! c2 < 0
  end do

  call Grid % Exchange_Cells_Real(turb % vis_t)
  call Grid % Exchange_Cells_Real(turb % vis_w)
  if(Flow % n_scalars > 0) call Grid % Exchange_Cells_Real(turb % diff_w)
  if(Flow % heat_transfer) call Grid % Exchange_Cells_Real(turb % con_w)

  end subroutine

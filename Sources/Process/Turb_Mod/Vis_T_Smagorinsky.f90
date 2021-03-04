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
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: u, v, w, t
  integer                   :: c, s, c1, c2
  real                      :: nx, ny, nz
  real                      :: cs, lf, u_tau_l, u_f, nc2
  real                      :: u_tan, a_pow, b_pow, nu, dely
  real                      :: beta, pr, ebf, u_plus, pr_t
!==============================================================================!

  ! Take aliases
  flow => turb % pnt_flow
  grid => flow % pnt_grid
  call Field_Mod_Alias_Momentum(flow, u, v, w)
  t    => flow % t

  !---------------!
  !               !
  !   SGS terms   !
  !               !
  !---------------!
  if(buoyancy) then
    call Field_Mod_Grad_Variable(flow, t)
  end if

  if(turb % model .eq. LES_SMAGORINSKY) then
    do c = 1, grid % n_cells
      lf = grid % vol(c)**ONE_THIRD

      ! if(nearest_wall_cell(c) .ne. 0) is needed for parallel version
      ! since the subdomains which do not "touch" wall
      ! has nearest_wall_cell(c) = 0. 
      if(turb % nearest_wall_cell(c) .ne. 0) then
        u_f = sqrt( flow % viscosity(c)                                &
                    * sqrt(  u % n(turb % nearest_wall_cell(c)) ** 2   &
                           + v % n(turb % nearest_wall_cell(c)) ** 2   &
                           + w % n(turb % nearest_wall_cell(c)) ** 2)  &
                   / (grid % wall_dist(turb % nearest_wall_cell(c))+TINY) )
        turb % y_plus(c) = grid % wall_dist(c) * u_f / flow % viscosity(c)
        cs = c_smag * (1.0 - exp(-turb % y_plus(c) / 25.0))
      else  
        cs = c_smag
      end if
      turb % vis_t(c) = flow % density(c)  &
                      * (lf*lf)            &  ! delta^2
                      * (cs*cs)            &  ! cs^2
                      * flow % shear(c)
    end do

  else if(turb % model .eq. LES_DYNAMIC) then
    do c = 1, grid % n_cells
      lf = grid % vol(c) ** ONE_THIRD
      turb % vis_t(c) = flow % density(c)  &
                      * (lf*lf)            &  ! delta^2
                      * turb % c_dyn(c)    &  ! c_dynamic
                      * flow % shear(c)
    end do
  else if(turb % model .eq. LES_WALE) then
    do c = 1, grid % n_cells
      lf = grid % vol(c)**ONE_THIRD
      turb % vis_t(c) = flow % density(c)  &
                      * (lf*lf)            &  ! delta^2
                      * (0.5*0.5)          &  ! cs^2
                      * turb % wale_v(c)
    end do
  end if

  if(buoyancy) then
    do c = 1, grid % n_cells
      nc2 = - flow % beta * (  grav_x * t % x(c)   &   
                             + grav_y * t % y(c)   &   
                             + grav_z * t % z(c))  
      turb % vis_t(c) = turb % vis_t(c)            &
             * max(1 - 2.5 * nc2 / (flow % shear(c) + TINY), 0.0)
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
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2  < 0) then

      nx = grid % sx(s) / grid % s(s)
      ny = grid % sy(s) / grid % s(s)
      nz = grid % sz(s) / grid % s(s)

      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

        u_tan = Field_Mod_U_Tan(flow, s)

        a_pow = 8.3
        b_pow = 1.0/7.0
        nu = flow % viscosity(c1) / flow % density(c1)
        dely = grid % wall_dist(c1)

        ! Calculate u_tau_l
        u_tau_l = ( u_tan/a_pow * (nu/dely)**b_pow ) ** (1.0/(1.0+b_pow))

        ! Calculate tau_wall (but it is never used)
        turb % tau_wall(c1) = flow % viscosity(c1) * u_tan / dely

        ! Calculate y+
        turb % y_plus(c1)  = dely * u_tau_l / nu
        if(turb % y_plus(c1)  >=  11.81) then
          ! This one is effective viscosity
          turb % vis_w(c1) = flow % density(c1) * u_tau_l * u_tau_l * dely  &
                           / abs(u_tan)
        else
          turb % vis_w(c1) = flow % viscosity(c1)                &
                        +      grid % fw(s)  * turb % vis_t(c1)  &
                        + (1.0-grid % fw(s)) * turb % vis_t(c2)
        end if

        if(flow % heat_transfer) then
          u_plus = u_tan / u_tau_l 
          pr_t = Turb_Mod_Prandtl_Number(turb, c1) 
          pr   = Field_Mod_Prandtl_Number(flow, c1)  ! laminar Prandtl number
          beta = 9.24 * ((pr/pr_t)**0.75 - 1.0)     &
               * (1.0 + 0.28 * exp(-0.007*pr/pr_t))
          ebf = Turb_Mod_Ebf_Scalar(turb, c1, pr) 
          turb % con_w(c1) =    turb % y_plus(c1)                         &
                              * flow % viscosity(c1)                      &
                              * flow % capacity(c1)                       &
                      / (  turb % y_plus(c1) * pr * exp(-1.0 * ebf)       &
                         + (u_plus + beta) * pr_t * exp(-1.0 / ebf) + TINY)
        end if

      end if  ! Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL or WALLFL
    end if    ! c2 < 0
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, turb % vis_t)
  call Grid_Mod_Exchange_Cells_Real(grid, turb % vis_w)

  end subroutine

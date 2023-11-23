!==============================================================================!
  subroutine Vis_T_Subgrid(Turb)
!------------------------------------------------------------------------------!
!   Calculates SGS stresses and turbulent viscosity for 'LES'.                 !
!------------------------------------------------------------------------------!
!   near(c) is the number of corresponding cell on the nearest wall.           !
!   In case that, in parallel executions, the subdomain does not have          !
!   any nearwall cells, the nearest_wall_cell(c) is zero.                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!------------------------------[Local parameters]------------------------------!
  real, parameter :: A_POW = 8.3
  real, parameter :: B_POW = 1.0/7.0
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w, t
  integer                   :: c, s, c1, c2
  real                      :: nx, ny, nz
  real                      :: cs, lf, u_tau, nc2, u_tan, nu
  real                      :: beta, pr, ebf, u_plus, pr_t, sc, z_o, kin_vis
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
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

  if(Turb % model .eq. LES_SMAGORINSKY) then
    do c = Cells_In_Domain_And_Buffers()
      lf = Grid % vol(c)**ONE_THIRD

      nu   = Flow % viscosity(c) / Flow % density(c)

      ! Tangential velocity.  Here it assumed that, as you approach the
      ! wall, the tangential velocity component is dominant and that the
      ! magnitude of velocity is close to tangential component.
      u_tan = sqrt(u % n(c)**2 + v % n(c)**2 + w % n(c)**2)

      ! Calculate u_tau, y+ and perform Van Driest damping
      u_tau = (u_tan/A_POW * (nu/Grid % wall_dist(c))**B_POW)   &
              ** (1.0/(1.0+B_POW))
      Turb % y_plus(c) = Grid % wall_dist(c) * u_tau / Flow % viscosity(c)
      cs = c_smag * (1.0 - exp(-Turb % y_plus(c) / 25.0))

      Turb % vis_t(c) = Flow % density(c)  &
                      * (lf*lf)            &  ! delta^2
                      * (cs*cs)            &  ! cs^2
                      * Flow % shear(c)
    end do

  else if(Turb % model .eq. LES_DYNAMIC) then
    do c = Cells_In_Domain_And_Buffers()
      lf = Grid % vol(c) ** ONE_THIRD
      Turb % vis_t(c) = Flow % density(c)  &
                      * (lf*lf)            &  ! delta^2
                      * Turb % c_dyn(c)    &  ! c_dynamic
                      * Flow % shear(c)
    end do

  else if(Turb % model .eq. LES_WALE) then
    do c = Cells_In_Domain_And_Buffers()
      lf = Grid % vol(c)**ONE_THIRD
      Turb % vis_t(c) = Flow % density(c)  &
                      * (lf*lf)            &  ! delta^2
                      * (0.5*0.5)          &  ! cs^2
                      * Turb % wale_v(c)
    end do
  end if

  !-------------------------------------------------!
  !   Modification of turbulent viscosity due to    !
  !   buoyancy according to Eidson, T., JFM, 1985   !
  !-------------------------------------------------!
  if(Flow % buoyancy .eq. THERMALLY_DRIVEN) then
    do c = Cells_In_Domain_And_Buffers()
      nc2 = max(- Flow % beta * (  Flow % grav_x * t % x(c)   &
                                 + Flow % grav_y * t % y(c)   &
                                 + Flow % grav_z * t % z(c)), 0.0)
      Turb % vis_t(c) = Turb % vis_t(c)            &
                      * sqrt(1 + 2.5 * nc2 / (Flow % shear(c)**2 + TINY))
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

        kin_vis =  Flow % viscosity(c1) / Flow % density(c1)

        ! Set up roughness coefficient
        z_o = Turb % Roughness_Coeff(c1, c2)

        u_tan = Flow % U_Tan(s)

        nu = Flow % viscosity(c1) / Flow % density(c1)

        ! Calculate u_tau for smooth wall
        u_tau = (u_tan/A_POW * (nu/Grid % wall_dist(c1))**B_POW) &
                ** (1.0/(1.0+B_POW))

        ! Calculate u_tau for rough wall
        if(z_o .gt. TINY) then
          u_tau = u_tan * kappa/log(Grid % wall_dist(c1)/z_o)
        end if

        ! Calculate u_tau according to Monin-Obukov Similarity Theory
        if(Flow % heat_transfer.and.Turb % monin_obukov) then
          u_tau = u_tau * Turb % Monin_Obukov_Momentum(abs(u_tan),  &
                  Grid % wall_dist(c1), z_o, t % n(c1),             &
                  t % n(c2), abs(Flow % grav_z))
        end if

        ! Calculate y+
        Turb % y_plus(c1) = Turb % Y_Plus_Rough_Walls(    &
                                   u_tau,                 &
                                   Grid % wall_dist(c1),  &
                                   kin_vis,               &
                                   z_o)

        u_plus = Turb % U_Plus_Log_Law(               &
                               Grid % wall_dist(c1),  &
                               Turb % y_plus(c1),     &
                               z_o)

        ! Effective viscosity above and below 11.18 threshold
        if(Turb % y_plus(c1)  >=  11.81) then
          Turb % vis_w(c1) = Flow % density(c1) * u_tau * u_tau  &
                             * Grid % wall_dist(c1) / abs(u_tan)
        else
          Turb % vis_w(c1) = Flow % viscosity(c1)                &
                        +      Grid % fw(s)  * Turb % vis_t(c1)  &
                        + (1.0-Grid % fw(s)) * Turb % vis_t(c2)
        end if

        if(Flow % heat_transfer) then
          pr_t = Turb % Prandtl_Turb(c1)
          pr   = Flow % Prandtl_Numb(c1)          ! laminar Prandtl number
          beta = Turb % Beta_Scalar(pr, pr_t)
          ! According to Toparlar et al. 2019 paper
          ! "CFD simulation of the near-neutral atmospheric boundary layer: New
          ! temperature inlet profile consistent with wall functions"
          if(z_o .gt. TINY) then
            beta = 0.0
          end if

          ebf = Turb % Ebf_Scalar(c1, pr)
          Turb % con_w(c1) =    Turb % y_plus(c1)                         &
                              * Flow % viscosity(c1)                      &
                              * Flow % capacity(c1)                       &
                      / (  Turb % y_plus(c1) * pr * exp(-1.0 * ebf)       &
                         + (u_plus + beta) * pr_t * exp(-1.0 / ebf) + TINY)

          if(Turb % monin_obukov) then
            Turb % con_w(c1) = pr_t * Turb % con_w(c1)                  &
                             * Turb % Monin_Obukov_Thermal(abs(u_tan),  &
                               Grid % wall_dist(c), z_o, t % n(c1),     &
                               t % n(c2),                               &
                               abs(Flow % grav_z))
          end if
        end if

        if(Flow % n_scalars > 0) then
          sc   = Flow % Schmidt_Numb(c1)          ! laminar Schmidt number
          beta = Turb % Beta_Scalar(sc, sc_t)
          ! According to Toparlar et al. 2019 paper
          ! "CFD simulation of the near-neutral atmospheric boundary layer: New
          ! temperature inlet profile consistent with wall functions"
          if(z_o .gt. TINY) then
            beta = 0.0
          end if

          ebf = Turb % Ebf_Scalar(c1, sc)
          Turb % diff_w(c1) =  Turb % y_plus(c1)                  &
              * (Flow % viscosity(c1)/Flow % density(c1))         &
              / (Turb % y_plus(c1) * sc * exp(-1.0 * ebf)         &
               + (u_plus + beta) * sc_t * exp(-1.0 / ebf) + TINY)

          if(Turb % monin_obukov) then
            Turb % diff_w(c1) = sc_t * Turb % diff_w(c1)                   &
                             * Turb % Monin_Obukov_Thermal(abs(u_tan),     &
                               Grid % wall_dist(c), z_o, t % n(c1), t % n(c2), &
                               abs(Flow % grav_z))
          end if
        end if

      end if  ! Grid % Bnd_Cond_Type(c2) .eq. WALL or WALLFL
    end if    ! c2 < 0
  end do

  call Grid % Exchange_Cells_Real(Turb % vis_t)
  call Grid % Exchange_Cells_Real(Turb % vis_w)
  if(Flow % n_scalars > 0) call Grid % Exchange_Cells_Real(Turb % diff_w)
  if(Flow % heat_transfer) call Grid % Exchange_Cells_Real(Turb % con_w)

  end subroutine

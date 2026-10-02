!==============================================================================!
  subroutine Wall_Function(Turb)
!------------------------------------------------------------------------------!
!   Computes turbulent viscosity for near-wall regions with wall functions     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w, t, kin
  type(Var_Type),   pointer :: vis
  integer                   :: s, c1, c2, reg
  real                      :: u_tau, u_tan, nu
  real                      :: beta, pr, ebf, u_plus, pr_t, sc, z_o, kin_vis
  real                      :: sc_tf
!------------------------------[Local parameters]------------------------------!
  real, parameter           :: A_POW = 8.3
  real, parameter           :: B_POW = 1.0/7.0
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  vis  => Turb % vis
  kin  => Turb % kin
  call Flow % Alias_Momentum(u, v, w)
  call Flow % Alias_Energy  (t)

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
  !---------------------------------------------!
  do reg = Boundary_Regions()

    if(Grid % region % type(reg) .eq. WALL .or.  &
       Grid % region % type(reg) .eq. WALLFL) then

      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        kin_vis =  Flow % viscosity(c1) / Flow % density(c1)

        ! Set up roughness coefficient
        z_o = Turb % Roughness_Coeff(c1, c2)

        u_tan = Flow % U_Tan(s)

        nu = Flow % viscosity(c1) / Flow % density(c1)

        ! Calculate u_tau for smooth wall
        if(Turb % model == K_EPS .or. Turb % model == K_EPS_ZETA_F .or. &
           Turb % model == HYBRID_LES_RANS .or. Turb % model == K_OMEGA_SST) then
          u_tau = Turb % c_mu25 * sqrt(kin % n(c1))
        else
          u_tau = Turb % U_Tau_Log_Law(u_tan,                &
                                       Grid % wall_dist(c1), &
                                       kin_vis,              &
                                       z_o)
        end if

        ! Calculate u_tau according to Monin-Obukov Similarity Theory
        if(Flow % heat_transfer .and. Turb % monin_obukov) then
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


        ebf = Turb % Ebf_Momentum(c1)

        if(Turb % y_plus(c1) < 3.0) then
          Turb % vis_w(c1) = Flow % viscosity(c1)
        else
          Turb % vis_w(c1) =    Turb % y_plus(c1) * Flow % viscosity(c1)  &
                           / (  Turb % y_plus(c1) * exp(-1.0 * ebf)       &
                              + u_plus * exp(-1.0/ebf) + TINY)
        end if

        if(Flow % heat_transfer) then
          pr_t = Turb % Prandtl_Turb(c1)
          pr   = Flow % Prandtl_Numb(c1)          ! laminar Prandtl number
          beta = Turb % Beta_Scalar(pr, pr_t)
          ! According to Toparlar et al. 2019 paper
          ! "CFD simulation of the near-neutral atmospheric boundary layer:
          ! New temperature inlet profile consistent with wall functions"
          if(z_o .gt. TINY) then
            beta = 0.0
          end if

          ! For y+ tending to zero, the expression below tends to molecular
          ! conductivity, but gives zero for y+ exactly zero.  With fluid at
          ! rest (y+ is zero at all walls) no heat could then enter through
          ! a wall with prescribed temperature, and a flow driven by buoyancy
          ! alone would never start.  Hence, the limit is set explicitly.
          if(Turb % y_plus(c1) < NANO) then
            Turb % con_w(c1) = Flow % conductivity(c1)
          else
            ebf = Turb % Ebf_Scalar(c1, pr)
            Turb % con_w(c1) =    Turb % y_plus(c1)                         &
                                * Flow % viscosity(c1)                      &
                                * Flow % capacity(c1)                       &
                        / (  Turb % y_plus(c1) * pr * exp(-1.0 * ebf)       &
                           + (u_plus + beta) * pr_t * exp(-1.0 / ebf) + TINY)
          end if

          if(Turb % monin_obukov) then
            Turb % con_w(c1) = pr_t * Turb % con_w(c1)                  &
                             * Turb % Monin_Obukov_Thermal(abs(u_tan),  &
                               Grid % wall_dist(c1), z_o, t % n(c1),    &
                               t % n(c2),                               &
                               abs(Flow % grav_z))
          end if
        end if  ! heat transfer

        if(Flow % n_scalars > 0) then
          sc   = Flow % Schmidt_Numb(c1)          ! laminar Schmidt number

          ! Turbulent Schmidt number; for LES models the one for
          ! sub-grid scale fluxes, as Prandtl_Turb does for pr_t above
          sc_tf = sc_t
          if(Turb % Les()) sc_tf = SC_T_LES

          beta = Turb % Beta_Scalar(sc, sc_tf)
          ! According to Toparlar et al. 2019 paper
          ! "CFD simulation of the near-neutral atmospheric boundary layer:
          ! New temperature inlet profile consistent with wall functions"
          if(z_o .gt. TINY) then
            beta = 0.0
          end if

          ! The same limit as for con_w above: molecular diffusivity
          if(Turb % y_plus(c1) < NANO) then
            Turb % diff_w(c1) = (Flow % viscosity(c1)/Flow % density(c1)) / sc
          else
            ebf = Turb % Ebf_Scalar(c1, sc)
            Turb % diff_w(c1) =  Turb % y_plus(c1)                  &
                * (Flow % viscosity(c1)/Flow % density(c1))         &
                / (Turb % y_plus(c1) * sc * exp(-1.0 * ebf)         &
                 + (u_plus + beta) * sc_tf * exp(-1.0 / ebf) + TINY)
          end if

          if(Turb % monin_obukov) then
            Turb % diff_w(c1) = sc_tf * Turb % diff_w(c1)                     &
                           * Turb % Monin_Obukov_Thermal(abs(u_tan),          &
                             Grid % wall_dist(c1), z_o, t % n(c1), t % n(c2), &
                             abs(Flow % grav_z))
          end if
        end if  ! n_scalars > 0

      end do
    end if      ! Grid % region % type(reg) .eq. WALL
  end do

  call Grid % Exchange_Cells_Real(Turb % vis_w)
  if(Flow % n_scalars > 0) then
    call Grid % Exchange_Cells_Real(Turb % diff_w)
  end if
  if(Flow % heat_transfer) then
    call Grid % Exchange_Cells_Real(Turb % con_w)
  end if

  end subroutine

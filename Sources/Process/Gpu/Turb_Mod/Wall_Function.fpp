!==============================================================================!
  subroutine Wall_Function(Turb, Grid, Flow)
!------------------------------------------------------------------------------!
!   Computes turbulent viscosity for near-wall regions with wall functions     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
  type(Grid_Type)          :: Grid
  type(Field_Type)         :: Flow
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type),   pointer :: u, v, w, t
  type(Var_Type),   pointer :: vis
  integer                   :: reg, s, c1, c2
  real                      :: u_tau, u_tan, nu
  real                      :: beta, pr, ebf, u_plus, pr_t, sc, z_o, kin_vis
!------------------------------[Local parameters]------------------------------!
  real, parameter           :: A_POW = 8.3
  real, parameter           :: B_POW = 1.0/7.0
!==============================================================================!

  ! Take aliases
  vis  => Turb % vis
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

    if(Grid % region % type(reg) .eq. WALL    .or.  &
       Grid % region % type(reg) .eq. WALLFL) then

      !$tf-acc loop begin
      do s = Faces_In_Region(reg)  ! all present
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        kin_vis =  Flow % viscosity(c1) / Flow % density(c1)

        ! Set up roughness coefficient
        z_o = Turb % Roughness_Coeff(Grid, c1, c2)

        u_tan = Flow % U_Tan(Grid, s)

        nu = Flow % viscosity(c1) / Flow % density(c1)

        ! Calculate u_tau for smooth wall
        u_tau = (u_tan/A_POW * (nu/Grid % wall_dist(c1))**B_POW) &
                ** (1.0/(1.0+B_POW))

        ! Calculate u_tau for rough wall
        if(z_o .gt. TINY) then
          u_tau = u_tan * Turb % kappa / log(Grid % wall_dist(c1)/z_o)
        end if

!1 !future ! Calculate u_tau according to Monin-Obukov Similarity Theory
!1 !future if(Flow % heat_transfer .and. Turb % monin_obukov) then
!1 !future   u_tau = u_tau * Turb % Monin_Obukov_Momentum(abs(u_tan),  &
!1 !future           Grid % wall_dist(c1), z_o, t % n(c1),             &
!1 !future           t % n(c2), abs(Flow % grav_z))
!1 !future end if

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
          Turb % vis_w(c1) = Turb % vis_t(c1) + Flow % viscosity(c1)
        else
          Turb % vis_w(c1) =    Turb % y_plus(c1) * Flow % viscosity(c1)  &
                           / (  Turb % y_plus(c1) * exp(-1.0 * ebf)       &
                              + u_plus * exp(-1.0/ebf) + TINY)
        end if

        if(Flow % heat_transfer) then
          pr_t = Turb % Prandtl_Turb(Flow, c1)
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

!1 !future   if(Turb % monin_obukov) then
!1 !future     Turb % con_w(c1) = pr_t * Turb % con_w(c1)                  &
!1 !future                      * Turb % Monin_Obukov_Thermal(abs(u_tan),  &
!1 !future                        Grid % wall_dist(c1), z_o, t % n(c1),    &
!1 !future                        t % n(c2),                               &
!1 !future                        abs(Flow % grav_z))
!1 !future   end if
        end if  ! heat transfer

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

!1 !future   if(Turb % monin_obukov) then
!1 !future     Turb % diff_w(c1) = sc_t * Turb % diff_w(c1)                      &
!1 !future                    * Turb % Monin_Obukov_Thermal(abs(u_tan),          &
!1 !future                      Grid % wall_dist(c1), z_o, t % n(c1), t % n(c2), &
!1 !future                      abs(Flow % grav_z))
!1 !future   end if
        end if  ! n_scalars > 0

      end do    ! s, faces
      !$tf-acc loop end

    end if      ! is it WALL or WALLFL?
  end do        ! through regions

  call Grid % Exchange_Cells_Real(Turb % vis_w)
  if(Flow % n_scalars > 0) call Grid % Exchange_Cells_Real(Turb % diff_w)
  if(Flow % heat_transfer) call Grid % Exchange_Cells_Real(Turb % con_w)

  end subroutine

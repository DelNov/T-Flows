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
  integer :: reg, s, c1, c2
  real    :: u_tau, u_tan, nu
  real    :: beta, pr, ebf, u_plus, pr_t, sc, z_o, kin_vis
!------------------------------[Local parameters]------------------------------!
  real, parameter :: A_POW = 8.3
  real, parameter :: B_POW = 1.0/7.0
!==============================================================================!

  !-------------------!
  !                   !
  !   Wall function   !
  !                   !
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

        !----------------------------------!
        !   Set up roughness coefficient   !
        !----------------------------------!
        z_o = Turb % z_o(c2)    ! take the value specified in control file
        if(z_o .gt. TINY) then  ! set lower limit based on wall distance
          z_o = max(Grid % wall_dist(c1)  &
              / (Turb % e_log * max(Turb % y_plus(c1), 1.0)), z_o)
        end if

        !------------------------------------------------!
        !   Work out the tangential velocity component   !
        !   (Assume it is the same ase the magnitude)    !
        !------------------------------------------------!
        u_tan = sqrt(  Flow % u % n(c1) ** 2   &
                     + Flow % v % n(c1) ** 2   &
                     + Flow % w % n(c1) ** 2)

        nu = Flow % viscosity(c1) / Flow % density(c1)

        !---------------------!
        !   Calculate u_tau   !
        !---------------------!

        ! First calculate u_tau for smooth wall ...
        u_tau = (u_tan/A_POW * (nu/Grid % wall_dist(c1))**B_POW) &
                ** (1.0/(1.0+B_POW))

        ! ... then correct, if needed, for rough wall
        if(z_o .gt. TINY) then
          u_tau = u_tan * Turb % kappa / log(Grid % wall_dist(c1)/z_o)
        end if

!1 !future ! Calculate u_tau according to Monin-Obukov Similarity Theory
!1 !future if(Flow % heat_transfer .and. Turb % monin_obukov) then
!1 !future   u_tau = u_tau * Turb % Monin_Obukov_Momentum(abs(u_tan),  &
!1 !future           Grid % wall_dist(c1), z_o, t % n(c1),             &
!1 !future           t % n(c2), abs(Flow % grav_z))
!1 !future end if

        !------------------!
        !   Calculate y+   !
        !------------------!
        Turb % y_plus(c1) = u_tau * (Grid % wall_dist(c1) + z_o) / kin_vis

        !------------------!
        !   Calculate u+   !
        !------------------!

        ! First salculate u+ for smooth walls ...
        u_plus = log( max(Turb % y_plus(c1), 1.05) * Turb % e_log )  &
               / Turb % kappa

        ! ... then correct, if needed, for rough walls
        if(z_o > TINY) then
          u_plus = log( (Grid % wall_dist(c1) + z_o) / z_o)  &
                 / (Turb % kappa + TINY) + TINY
        end if

        !-------------------------------------------------!
        !   Calculate blending coefficient for momentum   !
        !-------------------------------------------------!
        ebf = max(   0.01 * Turb % y_plus(c1)**4            &
                  / (1.0 + 5.0 * Turb % y_plus(c1)), TINY)

        if(Turb % y_plus(c1) < 3.0) then
          Turb % vis_w(c1) = Turb % vis_t(c1) + Flow % viscosity(c1)
        else
          Turb % vis_w(c1) =    Turb % y_plus(c1) * Flow % viscosity(c1)  &
                           / (  Turb % y_plus(c1) * exp(-1.0 * ebf)       &
                              + u_plus * exp(-1.0/ebf) + TINY)
        end if

        !-------------------------------------------------!
        !   For heat transfer problems, calculate con_w   !
        !-------------------------------------------------!
        if(Flow % heat_transfer) then
          pr_t = 0.4                      ! WARNING: GOOD ONLY FOR LES
          pr   = Flow % viscosity(c1)  &  ! laminar Prandtl number
               * Flow % capacity(c1)   &
               / Flow % conductivity(c1)
          beta = 9.24 * ((pr/pr_t)**0.75 - 1.0)  &
                      * (1.0 + 0.28 * exp(-0.007*pr/pr_t))
          ! According to Toparlar et al. 2019 paper
          ! "CFD simulation of the near-neutral atmospheric boundary layer: New
          ! temperature inlet profile consistent with wall functions"
          if(z_o .gt. TINY) then
            beta = 0.0
          end if

          ! Calculate blending coefficient for scalar(s)
          ebf = 0.01 * ((pr * Turb % y_plus(c1)) ** 4  &
              / ((1.0 + 5.0 * pr**3 * Turb % y_plus(c1)) + TINY))

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

        !-----------------------------------------------------!
        !   For scalar transport problems, calculate diff_w   !
        !-----------------------------------------------------!
        if(Flow % n_scalars > 0) then
          sc =  Flow % viscosity(c1)  &                       ! Schmidt number
             / (Flow % diffusivity(c1) * Flow % density(c1))
          beta = 9.24 * ((sc/sc_t)**0.75 - 1.0)  &
                      * (1.0 + 0.28 * exp(-0.007*sc/sc_t))
          ! According to Toparlar et al. 2019 paper
          ! "CFD simulation of the near-neutral atmospheric boundary layer: New
          ! temperature inlet profile consistent with wall functions"
          if(z_o .gt. TINY) then
            beta = 0.0
          end if

          ! Calculate blending coefficient for scalar(s)
          ebf = 0.01 * ((sc * Turb % y_plus(c1)) ** 4  &
              / ((1.0 + 5.0 * sc**3 * Turb % y_plus(c1)) + TINY))
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

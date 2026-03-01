!==============================================================================!
  subroutine Vis_T_Hybrid_Les_Prandtl(Turb)
!------------------------------------------------------------------------------!
!   Calculates SGS stresses and turbulent viscosity for 'LES'.                 !
!------------------------------------------------------------------------------!
!   near(c) is the number of corresponding cell on the nearest wall.
!   In case that, in parallel executions, the subdomain does not have
!   any nearwall cells, the near(c) is zero.
!   near(c) is calculated in NearWallCells.f90, only ones in the beginig
!   of a simulation.
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
  real                      :: fd    ! damping function
  real                      :: hwn   ! grid step size in wall-normal direction
  real                      :: hmax
  real                      :: cw, u_tan, u_tau, dw, lf_wm, kappa, dely, nu
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)
  t    => Flow % t

  !---------------!
  !   Constants   !  (Bad practice, constants should be in Turb_Mod.f90
  !---------------!
  cw    = 0.15 ! emperical constant [Shur 2008]
  kappa = 0.41

  !----------------------------!
  !   Model's eddy viscosity   !
  !----------------------------!
  do c = Cells_In_Domain_And_Buffers()

    hmax = Turb % h_max(c)
    hwn  = Turb % h_w(c)
    dw   = Grid % wall_dist(c)

    ! Wall-modeled LES length scale
    lf_wm = min(max(cw*dw,cw*hmax,hwn),hmax)

    nu   = Flow % viscosity(c) / Flow % density(c)
    dely = Grid % wall_dist(c)

    ! Tangential velocity.  Here it assumed that, as you approach the
    ! wall, the tangential velocity component is dominant and that the
    ! magnitude of velocity is close to tangential component.
    u_tan = sqrt(u % n(c)**2 + v % n(c)**2 + w % n(c)**2)

    ! Calculate u_tau, y+ and perform Van Driest damping
    u_tau = (u_tan/A_POW * (nu/dely)**B_POW) ** (1.0/(1.0+B_POW))
    Turb % y_plus(c) = Grid % wall_dist(c) * u_tau / Flow % viscosity(c)

    ! Piomelli damping function
    fd = 1.0 - exp(-(Turb % y_plus(c)/25.0)**3)

    ! Final SGS viscosity
    Turb % vis_t(c) = min((Turb % c_smag * lf_wm)**2, (kappa*dw)**2)  &
                    * Flow % shear(c) * fd

  end do

  !-----------------!
  !   Wall Region   !
  !-----------------+---------------------------!
  !  The procedure below calculates the vis..   !
  !  .. at the wall.                            !
  !----------------.----------------------------!
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    if(c2 < 0) then
      Turb % vis_w(c1) = Flow % viscosity(c1)            &
              +        Grid % fw(s)  * Turb % vis_t(c1)  &
              + (1.0 - Grid % fw(s)) * Turb % vis_t(c2)
    end if    ! c2 < 0
  end do

  call Grid % Exchange_Cells_Real(Turb % vis_t)
  call Grid % Exchange_Cells_Real(Turb % vis_w)

  end subroutine

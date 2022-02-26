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
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w, t
  integer                   :: c, s, c1, c2
  real                      :: fd    ! damping function
  real                      :: hwn   ! Grid step size in wall-normal direction 
  real                      :: hmax
  real                      :: cw, u_ff
  real                      :: dw
  real                      :: lf_wm
  real                      :: kappa
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

  ! Calculate model's eddy viscosity
  do c = 1, Grid % n_cells

    hmax = Turb % h_max(c)
    hwn  = Turb % h_w(c)
    dw   = Grid % wall_dist(c)

    ! Wall-modeled LES length scale
    lf_wm = min(max(cw*dw,cw*hmax,hwn),hmax)

    ! If(nearest_wall_cell(c) .ne. 0) is needed for parallel
    ! version since the subdomains which do not "touch" wall
    ! has nearest_wall_cell(c) = 0. 
    if(Turb % nearest_wall_cell(c) .ne. 0) then
      u_ff = sqrt( Flow % viscosity(c)  &
                  * sqrt(  u % n(Turb % nearest_wall_cell(c)) ** 2   &
                         + v % n(Turb % nearest_wall_cell(c)) ** 2   &
                         + w % n(Turb % nearest_wall_cell(c)) ** 2)  &
                 / (Grid % wall_dist(Turb % nearest_wall_cell(c))+TINY) )
      Turb % y_plus(c) = Grid % wall_dist(c) * u_ff / Flow % viscosity(c)

      ! Piomelli damping function
      fd = 1.0 - exp(-(Turb % y_plus(c)/25.0)**3)
    else
      fd = 1.0
    end if
    Turb % vis_t(c) = min((c_smag*lf_wm)**2, (kappa*dw)**2)  &
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

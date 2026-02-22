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
  integer                   :: c
  real                      :: cs, lf, u_tau, nc2, u_tan, nu
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
      cs = Turb % c_smag * (1.0 - exp(-Turb % y_plus(c) / 25.0))

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
      Turb % vis_t(c) = Flow % density(c)                &
                      * (lf*lf)                          &  ! delta^2
                      * (Turb % c_wale * Turb % c_wale)  &  ! cs^2
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
  !-------------------!
  call Turb % Wall_Function()

  call Grid % Exchange_Cells_Real(Turb % vis_t)

  end subroutine

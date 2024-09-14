!==============================================================================!
  subroutine Vis_T_Subgrid(Turb, Flow, Grid)
!------------------------------------------------------------------------------!
!>  Calculates SGS stresses and turbulent viscosity for 'LES'.                 !
!------------------------------------------------------------------------------!
!   near(c) is the number of corresponding cell on the nearest wall.           !
!   In case that, in parallel executions, the subdomain does not have          !
!   any nearwall cells, the nearest_wall_cell(c) is zero.                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb  !! parent class
  type(Field_Type), target :: Flow
  type(Grid_Type),  target :: Grid
!------------------------------[Local parameters]------------------------------!
  real, parameter :: A_POW = 8.3
  real, parameter :: B_POW = 1.0/7.0
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, s, c1, c2
  real    :: nx, ny, nz
  real    :: cs, lf, u_tau, nc2, u_tan, nu
  real    :: beta, pr, ebf, u_plus, pr_t, sc, z_o, kin_vis
!==============================================================================!

  !---------------!
  !               !
  !   SGS terms   !
  !               !
  !---------------!
  if(Turb % model .eq. LES_SMAGORINSKY) then

    !$acc parallel loop independent  &
    !$acc present(  &
    !$acc   grid_region_f_cell,  &
    !$acc   grid_region_l_cell,  &
    !$acc   grid_vol,  &
    !$acc   flow_viscosity,  &
    !$acc   flow_density,  &
    !$acc   flow_u_n,  &
    !$acc   flow_v_n,  &
    !$acc   flow_w_n,  &
    !$acc   grid_wall_dist,  &
    !$acc   turb_y_plus,  &
    !$acc   turb_vis_t,  &
    !$acc   flow_shear   &
    !$acc )
    do c = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions+1)
      lf = grid_vol(c)**ONE_THIRD

      nu = flow_viscosity(c) / flow_density(c)

      ! Tangential velocity.  Here it assumed that, as you approach the
      ! wall, the tangential velocity component is dominant and that the
      ! magnitude of velocity is close to tangential component.
      u_tan = sqrt(  flow_u_n(c)**2   &
                   + flow_v_n(c)**2   &
                   + flow_w_n(c)**2)

      ! Calculate u_tau, y+ and perform Van Driest damping
      u_tau = (u_tan/A_POW * (nu/grid_wall_dist(c))**B_POW)   &
              ** (1.0/(1.0+B_POW))
      turb_y_plus(c) = grid_wall_dist(c) * u_tau / flow_viscosity(c)
      cs = c_smag * (1.0 - exp(-turb_y_plus(c) / 25.0))

      turb_vis_t(c) = flow_density(c)  &
                      * (lf*lf)            &  ! delta^2
                      * (cs*cs)            &  ! cs^2
                      * flow_shear(c)
    end do
    !$acc end parallel

  else if(Turb % model .eq. LES_WALE) then

    !$acc parallel loop independent  &
    !$acc present(  &
    !$acc   grid_region_f_cell,  &
    !$acc   grid_region_l_cell,  &
    !$acc   grid_vol,  &
    !$acc   turb_vis_t,  &
    !$acc   flow_density,  &
    !$acc   turb_wale_v   &
    !$acc )
    do c = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions+1)
      lf = grid_vol(c)**ONE_THIRD
      turb_vis_t(c) = flow_density(c)  &
                      * (lf*lf)            &  ! delta^2
                      * (0.5*0.5)          &  ! cs^2
                      * turb_wale_v(c)
    end do
    !$acc end parallel

  end if

  call Grid % Exchange_Cells_Real(Turb % vis_t)

  end subroutine

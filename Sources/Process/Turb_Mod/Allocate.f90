!==============================================================================!
  subroutine Turb_Mod_Allocate(turb, flow)
!------------------------------------------------------------------------------!
!   Allocates memory for variables. It is called either from LoaRes            !
!   or from Processor.                                                         !
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Field_Mod,   only: Field_Type, heat_transfer, buoyancy
  use Grid_Mod,    only: Grid_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type),  target :: turb
  type(Field_Type), target :: flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: nb, nc
!==============================================================================!

  ! Store pointers
  turb % pnt_flow => flow
  turb % pnt_grid => flow % pnt_grid

  ! Take aliases
  grid => flow % pnt_grid
  nb = grid % n_bnd_cells
  nc = grid % n_cells

  ! Allocate deltas
  allocate(h_max(-nb:nc));  h_max = 0.
  allocate(h_min(-nb:nc));  h_min = 0.
  allocate(h_w  (-nb:nc));  h_w   = 0.

  !-----------------!
  !   K-eps model   !
  !-----------------!
  if(turbulence_model .eq. K_EPS) then

    ! Variables we solve for: k and epsilon
    call Var_Mod_Allocate_Solution('KIN', '', turb % kin, grid)
    call Var_Mod_Allocate_Solution('EPS', '', turb % eps, grid)

    ! Other turbulent quantities
    allocate(u_tau   (-nb:nc));  u_tau    = 0.
    allocate(p_kin   (-nb:nc));  p_kin    = 0.
    allocate(vis_t   (-nb:nc));  vis_t    = 0.
    allocate(vis_wall(-nb:nc));  vis_wall = 0.
    allocate(tau_wall(nc));                      tau_wall = 0.
    allocate(y_plus  (-nb:nc));  y_plus   = 0.

    if(heat_transfer) then
      call Var_Mod_Allocate_New_Only('UT', turb % ut, grid)
      call Var_Mod_Allocate_New_Only('VT', turb % vt, grid)
      call Var_Mod_Allocate_New_Only('WT', turb % wt, grid)
      allocate(con_wall(-nb:nc)); con_wall = 0.
    end if ! heat_transfer

    ! Turbulent statistics; if needed
    if(turbulence_statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(turb % u_mean(-nb:nc));  turb % u_mean = 0.
      allocate(turb % v_mean(-nb:nc));  turb % v_mean = 0.
      allocate(turb % w_mean(-nb:nc));  turb % w_mean = 0.
      allocate(turb % p_mean(-nb:nc));  turb % p_mean = 0.
      if(heat_transfer) then
        allocate(turb % t_mean(-nb:nc));  turb % t_mean = 0.
      end if

      ! Time-averaged modeled quantities
      allocate(turb % kin_mean(-nb:nc));  turb % kin_mean = 0.
      allocate(turb % eps_mean(-nb:nc));  turb % eps_mean = 0.

      ! Resolved Reynolds stresses
      allocate(turb % uu_res(-nb:nc));  turb % uu_res = 0.
      allocate(turb % vv_res(-nb:nc));  turb % vv_res = 0.
      allocate(turb % ww_res(-nb:nc));  turb % ww_res = 0.
      allocate(turb % uv_res(-nb:nc));  turb % uv_res = 0.
      allocate(turb % vw_res(-nb:nc));  turb % vw_res = 0.
      allocate(turb % uw_res(-nb:nc));  turb % uw_res = 0.

      ! Resolved turbulent heat fluxes
      if(heat_transfer) then
        allocate(turb % t2_res(-nb:nc));  turb % t2_res = 0.
        allocate(turb % ut_res(-nb:nc));  turb % ut_res = 0.
        allocate(turb % vt_res(-nb:nc));  turb % vt_res = 0.
        allocate(turb % wt_res(-nb:nc));  turb % wt_res = 0.
      end if ! heat_transfer

    end if ! turbulence_statistics

  end if ! K_EPS

  !------------------!
  !   K-eps-zeta-f   !
  !------------------!
  if(turbulence_model .eq. K_EPS_ZETA_F) then

    ! Main model's variables
    call Var_Mod_Allocate_Solution('KIN',  '', turb % kin,  grid)
    call Var_Mod_Allocate_Solution('EPS',  '', turb % eps,  grid)
    call Var_Mod_Allocate_Solution('ZETA', '', turb % zeta, grid)
    call Var_Mod_Allocate_Solution('F22',  '', turb % f22,  grid)

    ! Other variables such as time scale, length scale and production
    allocate(t_scale (-nb:nc));  t_scale  = 0.
    allocate(l_scale (-nb:nc));  l_scale  = 0.
    allocate(u_tau   (-nb:nc));  u_tau    = 0.
    allocate(p_kin   (-nb:nc));  p_kin    = 0.
    allocate(vis_t   (-nb:nc));  vis_t    = 0.
    allocate(vis_wall(-nb:nc));  vis_wall = 0.
    allocate(tau_wall(nc));                      tau_wall = 0.
    allocate(y_plus  (-nb:nc));  y_plus   = 0.

    ! Hydraulic roughness given by formula
    if(rough_walls) then
      allocate(z_o_f(-nb:nc));  z_o_f   = 0.
    end if

    if(heat_transfer) then
      call Var_Mod_Allocate_Solution('T2', '', turb % t2, grid)
      call Var_Mod_Allocate_New_Only('UT', turb % ut, grid)
      call Var_Mod_Allocate_New_Only('VT', turb % vt, grid)
      call Var_Mod_Allocate_New_Only('WT', turb % wt, grid)
      allocate(con_wall(-nb:nc)); con_wall = 0.
      allocate(p_t2    (-nb:nc)); p_t2     = 0.
    end if ! heat_transfer

    if(turbulence_statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(turb % u_mean(-nb:nc));  turb % u_mean = 0.
      allocate(turb % v_mean(-nb:nc));  turb % v_mean = 0.
      allocate(turb % w_mean(-nb:nc));  turb % w_mean = 0.
      allocate(turb % p_mean(-nb:nc));  turb % p_mean = 0.
      if(heat_transfer) then
        allocate(turb % t_mean(-nb:nc));  turb % t_mean = 0.
      end if

      ! Time-averaged modeled quantities
      allocate(turb % kin_mean (-nb:nc));  turb % kin_mean  = 0.
      allocate(turb % eps_mean (-nb:nc));  turb % eps_mean  = 0.
      allocate(turb % f22_mean (-nb:nc));  turb % f22_mean  = 0.
      allocate(turb % zeta_mean(-nb:nc));  turb % zeta_mean = 0.

      ! Resolved Reynolds stresses
      allocate(turb % uu_res(-nb:nc));  turb % uu_res = 0.
      allocate(turb % vv_res(-nb:nc));  turb % vv_res = 0.
      allocate(turb % ww_res(-nb:nc));  turb % ww_res = 0.
      allocate(turb % uv_res(-nb:nc));  turb % uv_res = 0.
      allocate(turb % vw_res(-nb:nc));  turb % vw_res = 0.
      allocate(turb % uw_res(-nb:nc));  turb % uw_res = 0.

      ! Resolved turbulent heat fluxes
      if(heat_transfer) then
        allocate(turb % t2_res(-nb:nc));  turb % t2_res = 0.
        allocate(turb % ut_res(-nb:nc));  turb % ut_res = 0.
        allocate(turb % vt_res(-nb:nc));  turb % vt_res = 0.
        allocate(turb % wt_res(-nb:nc));  turb % wt_res = 0.
      end if ! heat_transfer

    end if ! turbulence_statistics

    if(buoyancy) then
     allocate(g_buoy   (-nb:nc));  g_buoy    = 0.
     allocate(buoy_beta(-nb:nc));  buoy_beta = 0.
     allocate(g_kin    (-nb:nc));  g_kin     = 0.
    end if ! buoyancy

  end if ! K_EPS_ZETA_F

  !----------------------------!
  !   Reynolds stress models   !
  !----------------------------!
  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then

    ! Other variables such as time scale, length scale and production
    allocate(t_scale (-nb:nc));  t_scale  = 0.
    allocate(l_scale (-nb:nc));  l_scale  = 0.
    allocate(p_kin   (-nb:nc));  p_kin    = 0.
    allocate(u_tau   (-nb:nc));  u_tau    = 0.
    allocate(vis_t   (-nb:nc));  vis_t    = 0.
    allocate(vis_wall(-nb:nc));  vis_wall = 0.
    allocate(tau_wall(nc));                      tau_wall = 0.
    allocate(y_plus  (-nb:nc));  y_plus   = 0.

    ! Reynolds stresses
    call Var_Mod_Allocate_Solution('UU', '', turb % uu, grid)
    call Var_Mod_Allocate_Solution('VV', '', turb % vv, grid)
    call Var_Mod_Allocate_Solution('WW', '', turb % ww, grid)
    call Var_Mod_Allocate_Solution('UV', '', turb % uv, grid)
    call Var_Mod_Allocate_Solution('UW', '', turb % uw, grid)
    call Var_Mod_Allocate_Solution('VW', '', turb % vw, grid)

    call Var_Mod_Allocate_New_Only('KIN',     turb % kin, grid)
    call Var_Mod_Allocate_Solution('EPS', '', turb % eps, grid)

    if(heat_transfer) then
      call Var_Mod_Allocate_New_Only('UT', turb % ut, grid)
      call Var_Mod_Allocate_New_Only('VT', turb % vt, grid)
      call Var_Mod_Allocate_New_Only('WT', turb % wt, grid)
      allocate(con_wall(-nb:nc)); con_wall = 0.
    end if ! heat_transfer

    if(turbulence_statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(turb % u_mean(-nb:nc));  turb % u_mean = 0.
      allocate(turb % v_mean(-nb:nc));  turb % v_mean = 0.
      allocate(turb % w_mean(-nb:nc));  turb % w_mean = 0.
      allocate(turb % p_mean(-nb:nc));  turb % p_mean = 0.
      if(heat_transfer) then
        allocate(turb % t_mean(-nb:nc));  turb % t_mean = 0.
      end if

      ! Time-averaged modeled quantities
      allocate(turb % uu_mean (-nb:nc));  turb % uu_mean   = 0.
      allocate(turb % vv_mean (-nb:nc));  turb % vv_mean   = 0.
      allocate(turb % ww_mean (-nb:nc));  turb % ww_mean   = 0.
      allocate(turb % uv_mean (-nb:nc));  turb % uv_mean   = 0.
      allocate(turb % vw_mean (-nb:nc));  turb % vw_mean   = 0.
      allocate(turb % uw_mean (-nb:nc));  turb % uw_mean   = 0.
      allocate(turb % kin_mean(-nb:nc));  turb % kin_mean  = 0.
      allocate(turb % eps_mean(-nb:nc));  turb % eps_mean  = 0.

      ! Resolved Reynolds stresses
      allocate(turb % uu_res(-nb:nc));  turb % uu_res = 0.
      allocate(turb % vv_res(-nb:nc));  turb % vv_res = 0.
      allocate(turb % ww_res(-nb:nc));  turb % ww_res = 0.
      allocate(turb % uv_res(-nb:nc));  turb % uv_res = 0.
      allocate(turb % vw_res(-nb:nc));  turb % vw_res = 0.
      allocate(turb % uw_res(-nb:nc));  turb % uw_res = 0.

      ! Resolved turbulent heat fluxes
      if(heat_transfer) then
        allocate(turb % t2_res(-nb:nc));  turb % t2_res = 0.
        allocate(turb % ut_res(-nb:nc));  turb % ut_res = 0.
        allocate(turb % vt_res(-nb:nc));  turb % vt_res = 0.
        allocate(turb % wt_res(-nb:nc));  turb % wt_res = 0.
      end if ! heat_transfer

    end if ! turbulence_statistics

    if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then

      call Var_Mod_Allocate_Solution('F22', '', turb % f22, grid)

      if(turbulence_statistics) then
        allocate(turb % f22_mean(-nb:nc));  turb % f22_mean  = 0.
      end if ! turbulence_statistics

    end if ! RSM_MANCEAU_HANJALIC

  end if ! RSM_MANCEAU_HANJALIC & RSM_HANJALIC_JAKIRLIC

  !----------------------!
  !   Spalart Allmaras   !
  !----------------------!
  if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
     turbulence_model .eq. DES_SPALART) then

    call Var_Mod_Allocate_Solution('VIS', '', turb % vis, grid)

    ! Other variables such as time scale, length scale and production
    allocate(t_scale (-nb:nc));  t_scale  = 0.
    allocate(l_scale (-nb:nc));  l_scale  = 0.
    allocate(u_tau   (-nb:nc));  u_tau    = 0.
    allocate(p_kin   (-nb:nc));  p_kin    = 0.
    allocate(vis_t   (-nb:nc));  vis_t    = 0.
    allocate(vis_wall(-nb:nc));  vis_wall = 0.
    allocate(tau_wall(nc));                      tau_wall = 0.
    allocate(y_plus  (-nb:nc));  y_plus   = 0.

    if(heat_transfer) then
      call Var_Mod_Allocate_New_Only('UT', turb % ut, grid)
      call Var_Mod_Allocate_New_Only('VT', turb % vt, grid)
      call Var_Mod_Allocate_New_Only('WT', turb % wt, grid)
      allocate(con_wall(-nb:nc)); con_wall = 0.
    end if ! heat_transfer

    ! Turbulence statistics, if needed
    if(turbulence_statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(turb % u_mean(-nb:nc));  turb % u_mean = 0.
      allocate(turb % v_mean(-nb:nc));  turb % v_mean = 0.
      allocate(turb % w_mean(-nb:nc));  turb % w_mean = 0.
      allocate(turb % p_mean(-nb:nc));  turb % p_mean = 0.
      if(heat_transfer) then
        allocate(turb % t_mean(-nb:nc));  turb % t_mean = 0.
      end if

      ! Time-averaged modeled quantities
      allocate(turb % vis_mean (-nb:nc));  turb % vis_mean = 0.

      ! Resolved Reynolds stresses
      allocate(turb % uu_res(-nb:nc));  turb % uu_res = 0.
      allocate(turb % vv_res(-nb:nc));  turb % vv_res = 0.
      allocate(turb % ww_res(-nb:nc));  turb % ww_res = 0.
      allocate(turb % uv_res(-nb:nc));  turb % uv_res = 0.
      allocate(turb % vw_res(-nb:nc));  turb % vw_res = 0.
      allocate(turb % uw_res(-nb:nc));  turb % uw_res = 0.

      ! Resolved turbulent heat fluxes
      if(heat_transfer) then
        allocate(turb % t2_res(-nb:nc));  turb % t2_res = 0.
        allocate(turb % ut_res(-nb:nc));  turb % ut_res = 0.
        allocate(turb % vt_res(-nb:nc));  turb % vt_res = 0.
        allocate(turb % wt_res(-nb:nc));  turb % wt_res = 0.
      end if ! heat_transfer

    end if ! turbulence_statistics

  end if ! SPALART_ALLMARAS & DES_SPALART

  !-----------------------!
  !   Smagorinsky model   !
  !-----------------------!
  if(turbulence_model .eq. LES_SMAGORINSKY) then

    allocate(nearest_wall_cell(-nb:nc))
    nearest_wall_cell = 0

    ! Other variables such as time scale, length scale and production
    allocate(t_scale (-nb:nc));  t_scale  = 0.
    allocate(l_scale (-nb:nc));  l_scale  = 0.
    allocate(u_tau   (-nb:nc));  u_tau    = 0.
    allocate(p_kin   (-nb:nc));  p_kin    = 0.
    allocate(vis_t   (-nb:nc));  vis_t    = 0.
    allocate(vis_wall(-nb:nc));  vis_wall = 0.
    allocate(tau_wall(nc));                      tau_wall = 0.
    allocate(y_plus  (-nb:nc));  y_plus   = 0.

    if(heat_transfer) then
      call Var_Mod_Allocate_New_Only('UT', turb % ut, grid)
      call Var_Mod_Allocate_New_Only('VT', turb % vt, grid)
      call Var_Mod_Allocate_New_Only('WT', turb % wt, grid)
      allocate(con_wall(-nb:nc)); con_wall = 0.
    end if ! heat_transfer

    if(turbulence_statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(turb % u_mean(-nb:nc));  turb % u_mean = 0.
      allocate(turb % v_mean(-nb:nc));  turb % v_mean = 0.
      allocate(turb % w_mean(-nb:nc));  turb % w_mean = 0.
      allocate(turb % p_mean(-nb:nc));  turb % p_mean = 0.
      if(heat_transfer) then
        allocate(turb % t_mean(-nb:nc));  turb % t_mean = 0.
      end if

      ! Resolved Reynolds stresses
      allocate(turb % uu_res(-nb:nc));  turb % uu_res = 0.
      allocate(turb % vv_res(-nb:nc));  turb % vv_res = 0.
      allocate(turb % ww_res(-nb:nc));  turb % ww_res = 0.
      allocate(turb % uv_res(-nb:nc));  turb % uv_res = 0.
      allocate(turb % vw_res(-nb:nc));  turb % vw_res = 0.
      allocate(turb % uw_res(-nb:nc));  turb % uw_res = 0.

      ! Resolved turbulent heat fluxes
      if(heat_transfer) then
        allocate(turb % t2_res(-nb:nc));  turb % t2_res = 0.
        allocate(turb % ut_res(-nb:nc));  turb % ut_res = 0.
        allocate(turb % vt_res(-nb:nc));  turb % vt_res = 0.
        allocate(turb % wt_res(-nb:nc));  turb % wt_res = 0.
      end if ! heat_transfer

    end if ! turbulence_statistics

  end if ! LES_SMAGORINSKY

  !----------------!
  !   Wale model   !
  !----------------!
  if(turbulence_model .eq. LES_WALE) then

    allocate(wale_v(-nb:nc));  wale_v = 0.

    ! Other variables such as time scale, length scale and production
    allocate(t_scale (-nb:nc));  t_scale  = 0.
    allocate(l_scale (-nb:nc));  l_scale  = 0.
    allocate(u_tau   (-nb:nc));  u_tau    = 0.
    allocate(p_kin   (-nb:nc));  p_kin    = 0.
    allocate(vis_t   (-nb:nc));  vis_t    = 0.
    allocate(vis_wall(-nb:nc));  vis_wall = 0.
    allocate(tau_wall(nc));                      tau_wall = 0.
    allocate(y_plus  (-nb:nc));  y_plus   = 0.

    if(heat_transfer) then
      call Var_Mod_Allocate_New_Only('UT', turb % ut, grid)
      call Var_Mod_Allocate_New_Only('VT', turb % vt, grid)
      call Var_Mod_Allocate_New_Only('WT', turb % wt, grid)
      allocate(con_wall(-nb:nc)); con_wall = 0.
    end if ! heat_transfer

    if(turbulence_statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(turb % u_mean(-nb:nc));  turb % u_mean = 0.
      allocate(turb % v_mean(-nb:nc));  turb % v_mean = 0.
      allocate(turb % w_mean(-nb:nc));  turb % w_mean = 0.
      allocate(turb % p_mean(-nb:nc));  turb % p_mean = 0.
      if(heat_transfer) then
        allocate(turb % t_mean(-nb:nc));  turb % t_mean = 0.
      end if

      ! Resolved Reynolds stresses
      allocate(turb % uu_res(-nb:nc));  turb % uu_res = 0.
      allocate(turb % vv_res(-nb:nc));  turb % vv_res = 0.
      allocate(turb % ww_res(-nb:nc));  turb % ww_res = 0.
      allocate(turb % uv_res(-nb:nc));  turb % uv_res = 0.
      allocate(turb % vw_res(-nb:nc));  turb % vw_res = 0.
      allocate(turb % uw_res(-nb:nc));  turb % uw_res = 0.

      ! Resolved turbulent heat fluxes
      if(heat_transfer) then
        allocate(turb % t2_res(-nb:nc));  turb % t2_res = 0.
        allocate(turb % ut_res(-nb:nc));  turb % ut_res = 0.
        allocate(turb % vt_res(-nb:nc));  turb % vt_res = 0.
        allocate(turb % wt_res(-nb:nc));  turb % wt_res = 0.
      end if ! heat_transfer

    end if ! turbulence_statistics

  end if ! LES_WALE

  !-------------------!
  !   Dynamic model   !
  !-------------------!
  if(turbulence_model .eq. LES_DYNAMIC) then

    allocate(c_dyn(-nb:nc));  c_dyn = 0.

    ! Other variables such as time scale, length scale and production
    allocate(t_scale (-nb:nc));  t_scale  = 0.
    allocate(l_scale (-nb:nc));  l_scale  = 0.
    allocate(u_tau   (-nb:nc));  u_tau    = 0.
    allocate(p_kin   (-nb:nc));  p_kin    = 0.
    allocate(vis_t   (-nb:nc));  vis_t    = 0.
    allocate(vis_wall(-nb:nc));  vis_wall = 0.
    allocate(tau_wall(nc));                      tau_wall = 0.
    allocate(y_plus  (-nb:nc));  y_plus   = 0.

    if(heat_transfer) then

      call Var_Mod_Allocate_New_Only('UT', turb % ut, grid)
      call Var_Mod_Allocate_New_Only('VT', turb % vt, grid)
      call Var_Mod_Allocate_New_Only('WT', turb % wt, grid)
      allocate(con_wall(-nb:nc)); con_wall = 0.

    end if ! heat_transfer

    if(turbulence_statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(turb % u_mean(-nb:nc));  turb % u_mean = 0.
      allocate(turb % v_mean(-nb:nc));  turb % v_mean = 0.
      allocate(turb % w_mean(-nb:nc));  turb % w_mean = 0.
      allocate(turb % p_mean(-nb:nc));  turb % p_mean = 0.
      if(heat_transfer) then
        allocate(turb % t_mean(-nb:nc));  turb % t_mean = 0.
      end if

      ! Resolved Reynolds stresses
      allocate(turb % uu_res(-nb:nc));  turb % uu_res = 0.
      allocate(turb % vv_res(-nb:nc));  turb % vv_res = 0.
      allocate(turb % ww_res(-nb:nc));  turb % ww_res = 0.
      allocate(turb % uv_res(-nb:nc));  turb % uv_res = 0.
      allocate(turb % vw_res(-nb:nc));  turb % vw_res = 0.
      allocate(turb % uw_res(-nb:nc));  turb % uw_res = 0.

      ! Resolved turbulent heat fluxes
      if(heat_transfer) then
        allocate(turb % t2_res(-nb:nc));  turb % t2_res = 0.
        allocate(turb % ut_res(-nb:nc));  turb % ut_res = 0.
        allocate(turb % vt_res(-nb:nc));  turb % vt_res = 0.
        allocate(turb % wt_res(-nb:nc));  turb % wt_res = 0.
      end if ! heat_transfer

    end if ! turbulence_statistics

  end if ! LES_DYNAMIC

  !-------------------!
  !   Hybrid model    !
  !-------------------!
  if(turbulence_model .eq. HYBRID_LES_RANS) then

    ! Main model's variables (for RANS part)
    call Var_Mod_Allocate_Solution('KIN',  '', turb % kin,  grid)
    call Var_Mod_Allocate_Solution('EPS',  '', turb % eps,  grid)
    call Var_Mod_Allocate_Solution('ZETA', '', turb % zeta, grid)
    call Var_Mod_Allocate_Solution('F22',  '', turb % f22,  grid)

    ! Main model's variables (for LES part)
    allocate(c_dyn(-nb:nc));  c_dyn = 0.

    ! Other variables such as time scale, length scale and production
    allocate(vis_t_eff(-nb:nc));  vis_t_eff = 0.
    allocate(vis_t_sgs(-nb:nc));  vis_t_sgs = 0.
    allocate(t_scale  (-nb:nc));  t_scale   = 0.
    allocate(l_scale  (-nb:nc));  l_scale   = 0.
    allocate(u_tau    (-nb:nc));  u_tau     = 0.
    allocate(p_kin    (-nb:nc));  p_kin     = 0.
    allocate(vis_t    (-nb:nc));  vis_t     = 0.
    allocate(vis_wall (-nb:nc));  vis_wall  = 0.
    allocate(tau_wall (nc));                      tau_wall  = 0.
    allocate(y_plus   (-nb:nc));  y_plus    = 0.

    if(heat_transfer) then
      call Var_Mod_Allocate_Solution('T2', '', turb % t2, grid)
      call Var_Mod_Allocate_New_Only('UT', turb % ut, grid)
      call Var_Mod_Allocate_New_Only('VT', turb % vt, grid)
      call Var_Mod_Allocate_New_Only('WT', turb % wt, grid)
      allocate(con_wall(-nb:nc)); con_wall = 0.
    end if ! heat_transfer

    if(turbulence_statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(turb % u_mean(-nb:nc));  turb % u_mean = 0.
      allocate(turb % v_mean(-nb:nc));  turb % v_mean = 0.
      allocate(turb % w_mean(-nb:nc));  turb % w_mean = 0.
      allocate(turb % p_mean(-nb:nc));  turb % p_mean = 0.
      if(heat_transfer) then
        allocate(turb % t_mean(-nb:nc));  turb % t_mean = 0.
      end if

      ! Time-averaged modeled quantities
      allocate(turb % kin_mean (-nb:nc));  turb % kin_mean  = 0.
      allocate(turb % eps_mean (-nb:nc));  turb % eps_mean  = 0.
      allocate(turb % f22_mean (-nb:nc));  turb % f22_mean  = 0.
      allocate(turb % zeta_mean(-nb:nc));  turb % zeta_mean = 0.

      ! Resolved Reynolds stresses
      allocate(turb % uu_res(-nb:nc));  turb % uu_res = 0.
      allocate(turb % vv_res(-nb:nc));  turb % vv_res = 0.
      allocate(turb % ww_res(-nb:nc));  turb % ww_res = 0.
      allocate(turb % uv_res(-nb:nc));  turb % uv_res = 0.
      allocate(turb % vw_res(-nb:nc));  turb % vw_res = 0.
      allocate(turb % uw_res(-nb:nc));  turb % uw_res = 0.

      ! Resolved turbulent heat fluxes
      if(heat_transfer) then
        allocate(turb % t2_res(-nb:nc));  turb % t2_res = 0.
        allocate(turb % ut_res(-nb:nc));  turb % ut_res = 0.
        allocate(turb % vt_res(-nb:nc));  turb % vt_res = 0.
        allocate(turb % wt_res(-nb:nc));  turb % wt_res = 0.
      end if ! heat_transfer

    end if ! turbulence_statistics

    if(buoyancy) then
      allocate(g_buoy   (-nb:nc));  g_buoy    = 0.
      allocate(buoy_beta(-nb:nc));  buoy_beta = 0.
      allocate(g_kin    (-nb:nc));  g_kin     = 0.
    end if ! buoyancy

  end if ! HYBRID_LES_RANS

  end subroutine

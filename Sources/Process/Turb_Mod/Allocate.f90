!==============================================================================!
  subroutine Turb_Mod_Allocate(turb, Flow)
!------------------------------------------------------------------------------!
!   Allocates memory for variables. It is called either from LoaRes            !
!   or from Processor.                                                         !
!------------------------------------------------------------------------------!
!  use Multiphase_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type),  target :: turb
  type(Field_Type), target :: Flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: nb, nc
!==============================================================================!

  ! Store pointers
  turb % pnt_flow => Flow
  turb % pnt_grid => Flow % pnt_grid

  ! Take aliases
  Grid => Flow % pnt_grid
  nb = Grid % n_bnd_cells
  nc = Grid % n_cells

  ! Allocate deltas
  allocate(turb % h_max(-nb:nc));  turb % h_max = 0.
  allocate(turb % h_min(-nb:nc));  turb % h_min = 0.
  allocate(turb % h_w  (-nb:nc));  turb % h_w   = 0.

  allocate(turb % tau_wall(-nb:nc));  turb % tau_wall = 0.

  !-----------------!
  !   K-eps model   !
  !-----------------!
  if(turb % model .eq. K_EPS) then

    ! Variables we solve for: k and epsilon
    call Var_Mod_Allocate_Solution(turb % kin, Grid, 'KIN', '')
    call Var_Mod_Allocate_Solution(turb % eps, Grid, 'EPS', '')

    ! Other turbulent quantities
    allocate(turb % vis_t (-nb:nc));  turb % vis_t  = 0.
    allocate(turb % vis_w (-nb:nc));  turb % vis_w  = 0.  ! wall visc
    allocate(turb % p_kin (-nb:nc));  turb % p_kin  = 0.
    allocate(turb % y_plus(-nb:nc));  turb % y_plus = 0.

    if(Flow % heat_transfer) then
      call Var_Mod_Allocate_Solution(turb % t2, Grid, 'T2', '')
      call Var_Mod_Allocate_New_Only(turb % ut, Grid, 'UT')
      call Var_Mod_Allocate_New_Only(turb % vt, Grid, 'VT')
      call Var_Mod_Allocate_New_Only(turb % wt, Grid, 'WT')
      allocate(turb % con_w(-nb:nc));  turb % con_w = 0.  ! wall cond
      allocate(turb % p_t2 (-nb:nc));  turb % p_t2  = 0.

      ! Reynolds stresses
      call Var_Mod_Allocate_New_Only(turb % uu, Grid, 'UU')
      call Var_Mod_Allocate_New_Only(turb % vv, Grid, 'VV')
      call Var_Mod_Allocate_New_Only(turb % ww, Grid, 'WW')
      call Var_Mod_Allocate_New_Only(turb % uv, Grid, 'UV')
      call Var_Mod_Allocate_New_Only(turb % uw, Grid, 'UW')
      call Var_Mod_Allocate_New_Only(turb % vw, Grid, 'VW')
    end if ! Flow % heat_transfer

    !  Wall difussivity for user scalar
    if(Flow % n_scalars > 0) then            
      allocate(turb % diff_w(-nb:nc));  turb % diff_w = 0.  
    end if

    ! Turbulent statistics; if needed
    if(turb % statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(turb % u_mean(-nb:nc));  turb % u_mean = 0.
      allocate(turb % v_mean(-nb:nc));  turb % v_mean = 0.
      allocate(turb % w_mean(-nb:nc));  turb % w_mean = 0.
      allocate(turb % p_mean(-nb:nc));  turb % p_mean = 0.
      if(Flow % heat_transfer) then
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
      if(Flow % heat_transfer) then
        allocate(turb % t2_res (-nb:nc));  turb % t2_res  = 0.
        allocate(turb % ut_res (-nb:nc));  turb % ut_res  = 0.
        allocate(turb % vt_res (-nb:nc));  turb % vt_res  = 0.
        allocate(turb % wt_res (-nb:nc));  turb % wt_res  = 0.
        allocate(turb % ut_mean(-nb:nc));  turb % ut_mean = 0.
        allocate(turb % vt_mean(-nb:nc));  turb % vt_mean = 0.
        allocate(turb % wt_mean(-nb:nc));  turb % wt_mean = 0.
      end if ! Flow % heat_transfer

    end if ! turb % statistics

    if(Flow % buoyancy .eq. THERMALLY_DRIVEN) then
      allocate(turb % g_buoy(-nb:nc));  turb % g_buoy = 0.
    end if ! Flow % buoyancy .eq. THERMALLY_DRIVEN

  end if ! K_EPS

  !------------------!
  !   K-eps-zeta-f   !
  !------------------!
  if(turb % model .eq. K_EPS_ZETA_F) then

    ! Main model's variables
    call Var_Mod_Allocate_Solution(turb % kin,  Grid, 'KIN',  '')
    call Var_Mod_Allocate_Solution(turb % eps,  Grid, 'EPS',  '')
    call Var_Mod_Allocate_Solution(turb % zeta, Grid, 'ZETA', '')
    call Var_Mod_Allocate_Solution(turb % f22,  Grid, 'F22',  '')

    ! Other variables such as time scale, length scale and production
    allocate(turb % t_scale (-nb:nc));  turb % t_scale = 0.
    allocate(turb % l_scale (-nb:nc));  turb % l_scale = 0.
    allocate(turb % alpha_l (-nb:nc));  turb % alpha_l = 0.
    allocate(turb % alpha_u (-nb:nc));  turb % alpha_u = 0.
    allocate(turb % vis_t   (-nb:nc));  turb % vis_t   = 0.
    allocate(turb % vis_w   (-nb:nc));  turb % vis_w   = 0.  ! wall visc
    allocate(turb % p_kin   (-nb:nc));  turb % p_kin   = 0.
    allocate(turb % y_plus  (-nb:nc));  turb % y_plus  = 0.

    ! Hydraulic roughness given by formula
    if(turb % rough_walls) then
      allocate(turb % z_o_f  (-nb:nc)); turb % z_o_f   = -1.0
      allocate(turb % id_zone(-nb:nc)); turb % id_zone =  0.0
    end if

    ! Reynolds stresses
    call Var_Mod_Allocate_New_Only(turb % uu, Grid, 'UU')
    call Var_Mod_Allocate_New_Only(turb % vv, Grid, 'VV')
    call Var_Mod_Allocate_New_Only(turb % ww, Grid, 'WW')
    call Var_Mod_Allocate_New_Only(turb % uv, Grid, 'UV')
    call Var_Mod_Allocate_New_Only(turb % uw, Grid, 'UW')
    call Var_Mod_Allocate_New_Only(turb % vw, Grid, 'VW')

    if(Flow % heat_transfer) then
      call Var_Mod_Allocate_Solution(turb % t2, Grid, 'T2', '')
      call Var_Mod_Allocate_New_Only(turb % ut, Grid, 'UT')
      call Var_Mod_Allocate_New_Only(turb % vt, Grid, 'VT')
      call Var_Mod_Allocate_New_Only(turb % wt, Grid, 'WT')
      allocate(turb % con_w(-nb:nc));  turb % con_w = 0.  ! wall cond
      allocate(turb % p_t2 (-nb:nc));  turb % p_t2  = 0.

    end if ! Flow % heat_transfer

    !  Wall difussivity for user scalar
    if(Flow % n_scalars > 0) then            
      allocate(turb % diff_w(-nb:nc));  turb % diff_w = 0.  
    end if
   
    if(turb % statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(turb % u_mean(-nb:nc));  turb % u_mean = 0.
      allocate(turb % v_mean(-nb:nc));  turb % v_mean = 0.
      allocate(turb % w_mean(-nb:nc));  turb % w_mean = 0.
      allocate(turb % p_mean(-nb:nc));  turb % p_mean = 0.
      if(Flow % heat_transfer) then
        allocate(turb % t_mean(-nb:nc));  turb % t_mean = 0.
      end if

      ! Time-averaged modeled quantities
      allocate(turb % kin_mean (-nb:nc));  turb % kin_mean  = 0.
      allocate(turb % eps_mean (-nb:nc));  turb % eps_mean  = 0.
      allocate(turb % f22_mean (-nb:nc));  turb % f22_mean  = 0.
      allocate(turb % zeta_mean(-nb:nc));  turb % zeta_mean = 0.
      if(Flow % heat_transfer) then
        allocate(turb % t2_mean(-nb:nc));  turb % t2_mean   = 0.
      end if

      ! Resolved Reynolds stresses
      allocate(turb % uu_res(-nb:nc));  turb % uu_res = 0.
      allocate(turb % vv_res(-nb:nc));  turb % vv_res = 0.
      allocate(turb % ww_res(-nb:nc));  turb % ww_res = 0.
      allocate(turb % uv_res(-nb:nc));  turb % uv_res = 0.
      allocate(turb % vw_res(-nb:nc));  turb % vw_res = 0.
      allocate(turb % uw_res(-nb:nc));  turb % uw_res = 0.

      ! Resolved turbulent heat fluxes
      if(Flow % heat_transfer) then
        allocate(turb % t2_res (-nb:nc));  turb % t2_res  = 0.
        allocate(turb % ut_res (-nb:nc));  turb % ut_res  = 0.
        allocate(turb % vt_res (-nb:nc));  turb % vt_res  = 0.
        allocate(turb % wt_res (-nb:nc));  turb % wt_res  = 0.
        allocate(turb % ut_mean(-nb:nc));  turb % ut_mean = 0.
        allocate(turb % vt_mean(-nb:nc));  turb % vt_mean = 0.
        allocate(turb % wt_mean(-nb:nc));  turb % wt_mean = 0.
      end if ! Flow % heat_transfer

    end if ! turb % statistics

    if(Flow % buoyancy .eq. THERMALLY_DRIVEN) then
      allocate(turb % g_buoy(-nb:nc));  turb % g_buoy = 0.
    end if ! Flow % buoyancy .eq. THERMALLY_DRIVEN

  end if ! K_EPS_ZETA_F

  !----------------------------!
  !   Reynolds stress models   !
  !----------------------------!
  if(turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turb % model .eq. RSM_HANJALIC_JAKIRLIC) then

    ! Other variables such as time scale, length scale and production
    allocate(turb % t_scale(-nb:nc));  turb % t_scale = 0.
    allocate(turb % l_scale(-nb:nc));  turb % l_scale = 0.
    allocate(turb % vis_t  (-nb:nc));  turb % vis_t   = 0.
    allocate(turb % vis_w  (-nb:nc));  turb % vis_w   = 0.  ! wall visc
    allocate(turb % p_kin  (-nb:nc));  turb % p_kin   = 0.
    allocate(turb % y_plus (-nb:nc));  turb % y_plus  = 0.

    ! Reynolds stresses
    call Var_Mod_Allocate_Solution(turb % uu, Grid, 'UU', '')
    call Var_Mod_Allocate_Solution(turb % vv, Grid, 'VV', '')
    call Var_Mod_Allocate_Solution(turb % ww, Grid, 'WW', '')
    call Var_Mod_Allocate_Solution(turb % uv, Grid, 'UV', '')
    call Var_Mod_Allocate_Solution(turb % uw, Grid, 'UW', '')
    call Var_Mod_Allocate_Solution(turb % vw, Grid, 'VW', '')

    call Var_Mod_Allocate_New_Only(turb % kin, Grid, 'KIN')
    call Var_Mod_Allocate_Solution(turb % eps, Grid, 'EPS', '')

    if(Flow % heat_transfer) then
      call Var_Mod_Allocate_New_Only(turb % ut, Grid, 'UT')
      call Var_Mod_Allocate_New_Only(turb % vt, Grid, 'VT')
      call Var_Mod_Allocate_New_Only(turb % wt, Grid, 'WT')
      allocate(turb % con_w(-nb:nc));  turb % con_w = 0.  ! wall cond
    end if ! Flow % heat_transfer

    if(turb % statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(turb % u_mean(-nb:nc));  turb % u_mean = 0.
      allocate(turb % v_mean(-nb:nc));  turb % v_mean = 0.
      allocate(turb % w_mean(-nb:nc));  turb % w_mean = 0.
      allocate(turb % p_mean(-nb:nc));  turb % p_mean = 0.
      if(Flow % heat_transfer) then
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
      if(Flow % heat_transfer) then
        allocate(turb % t2_res (-nb:nc));  turb % t2_res  = 0.
        allocate(turb % ut_res (-nb:nc));  turb % ut_res  = 0.
        allocate(turb % vt_res (-nb:nc));  turb % vt_res  = 0.
        allocate(turb % wt_res (-nb:nc));  turb % wt_res  = 0.
        allocate(turb % ut_mean(-nb:nc));  turb % ut_mean = 0.
        allocate(turb % vt_mean(-nb:nc));  turb % vt_mean = 0.
        allocate(turb % wt_mean(-nb:nc));  turb % wt_mean = 0.
      end if ! Flow % heat_transfer

    end if ! turb % statistics

    if(turb % model .eq. RSM_MANCEAU_HANJALIC) then

      call Var_Mod_Allocate_Solution(turb % f22, Grid, 'F22', '')

      if(turb % statistics) then
        allocate(turb % f22_mean(-nb:nc));  turb % f22_mean  = 0.
      end if ! turb % statistics

    end if ! RSM_MANCEAU_HANJALIC

    if(Flow % buoyancy .eq. THERMALLY_DRIVEN) then
      allocate(turb % g_buoy(-nb:nc));  turb % g_buoy = 0.
    end if ! Flow % buoyancy .eq. THERMALLY_DRIVEN

  end if ! RSM_MANCEAU_HANJALIC & RSM_HANJALIC_JAKIRLIC

  !----------------------!
  !   Spalart Allmaras   !
  !----------------------!
  if(turb % model .eq. SPALART_ALLMARAS .or.  &
     turb % model .eq. DES_SPALART) then

    call Var_Mod_Allocate_Solution(turb % vis, Grid, 'VIS', '')

    ! Other variables such as time scale, length scale and production
    allocate(turb % t_scale (-nb:nc));  turb % t_scale = 0.
    allocate(turb % l_scale (-nb:nc));  turb % l_scale = 0.
    allocate(turb % vis_t   (-nb:nc));  turb % vis_t   = 0.
    allocate(turb % vis_w   (-nb:nc));  turb % vis_w   = 0.  ! wall visc
    allocate(turb % p_kin   (-nb:nc));  turb % p_kin   = 0.
    allocate(turb % y_plus  (-nb:nc));  turb % y_plus  = 0.

    if(Flow % heat_transfer) then
      call Var_Mod_Allocate_New_Only(turb % ut, Grid, 'UT')
      call Var_Mod_Allocate_New_Only(turb % vt, Grid, 'VT')
      call Var_Mod_Allocate_New_Only(turb % wt, Grid, 'WT')
      allocate(turb % con_w(-nb:nc));  turb % con_w = 0.  ! wall cond
    end if ! Flow % heat_transfer

    ! Turbulence statistics, if needed
    if(turb % statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(turb % u_mean(-nb:nc));  turb % u_mean = 0.
      allocate(turb % v_mean(-nb:nc));  turb % v_mean = 0.
      allocate(turb % w_mean(-nb:nc));  turb % w_mean = 0.
      allocate(turb % p_mean(-nb:nc));  turb % p_mean = 0.
      if(Flow % heat_transfer) then
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
      if(Flow % heat_transfer) then
        allocate(turb % t2_res(-nb:nc));  turb % t2_res = 0.
        allocate(turb % ut_res(-nb:nc));  turb % ut_res = 0.
        allocate(turb % vt_res(-nb:nc));  turb % vt_res = 0.
        allocate(turb % wt_res(-nb:nc));  turb % wt_res = 0.
      end if ! Flow % heat_transfer

    end if ! turb % statistics

  end if ! SPALART_ALLMARAS & DES_SPALART

  !-----------------------!
  !   Smagorinsky model   !
  !-----------------------!
  if(turb % model .eq. LES_SMAGORINSKY) then

    allocate(turb % nearest_wall_cell(-nb:nc))
    turb % nearest_wall_cell = 0

    ! Other variables such as time scale, length scale and production
    allocate(turb % t_scale (-nb:nc));  turb % t_scale = 0.
    allocate(turb % l_scale (-nb:nc));  turb % l_scale = 0.
    allocate(turb % vis_t   (-nb:nc));  turb % vis_t   = 0.
    allocate(turb % vis_w   (-nb:nc));  turb % vis_w   = 0.  ! wall visc
    allocate(turb % p_kin   (-nb:nc));  turb % p_kin   = 0.
    allocate(turb % y_plus  (-nb:nc));  turb % y_plus  = 0.

    if(Flow % heat_transfer) then
      call Var_Mod_Allocate_New_Only(turb % ut, Grid, 'UT')
      call Var_Mod_Allocate_New_Only(turb % vt, Grid, 'VT')
      call Var_Mod_Allocate_New_Only(turb % wt, Grid, 'WT')
      allocate(turb % con_w(-nb:nc));  turb % con_w = 0.  ! wall cond
    end if ! Flow % heat_transfer

    if(turb % statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(turb % u_mean(-nb:nc));  turb % u_mean = 0.
      allocate(turb % v_mean(-nb:nc));  turb % v_mean = 0.
      allocate(turb % w_mean(-nb:nc));  turb % w_mean = 0.
      allocate(turb % p_mean(-nb:nc));  turb % p_mean = 0.
      if(Flow % heat_transfer) then
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
      if(Flow % heat_transfer) then
        allocate(turb % t2_res(-nb:nc));  turb % t2_res = 0.
        allocate(turb % ut_res(-nb:nc));  turb % ut_res = 0.
        allocate(turb % vt_res(-nb:nc));  turb % vt_res = 0.
        allocate(turb % wt_res(-nb:nc));  turb % wt_res = 0.
      end if ! Flow % heat_transfer

    end if ! turb % statistics

  end if ! LES_SMAGORINSKY


  !------------------------------!
  !   Hybrid_Les_Prandtl model   !
  !------------------------------!
  if(turb % model .eq. HYBRID_LES_PRANDTL) then

    allocate(turb % nearest_wall_cell(-nb:nc))
    turb % nearest_wall_cell = 0
 
    ! Dynamic Smagorinsky constant for particle SGS models
    allocate(turb % c_dyn(-nb:nc));  turb % c_dyn = 0.

    ! Other variables such as time scale, length scale and production
    allocate(turb % t_scale (-nb:nc));  turb % t_scale = 0.
    allocate(turb % l_scale (-nb:nc));  turb % l_scale = 0.
    allocate(turb % vis_t   (-nb:nc));  turb % vis_t   = 0.
    allocate(turb % vis_w   (-nb:nc));  turb % vis_w   = 0.  ! wall visc
    allocate(turb % p_kin   (-nb:nc));  turb % p_kin   = 0.
    allocate(turb % y_plus  (-nb:nc));  turb % y_plus  = 0.

    if(Flow % heat_transfer) then
      call Var_Mod_Allocate_New_Only(turb % ut, Grid, 'UT')
      call Var_Mod_Allocate_New_Only(turb % vt, Grid, 'VT')
      call Var_Mod_Allocate_New_Only(turb % wt, Grid, 'WT')
      allocate(turb % con_w(-nb:nc));  turb % con_w = 0.  ! wall cond
    end if ! Flow % heat_transfer

    if(turb % statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(turb % u_mean(-nb:nc));  turb % u_mean = 0.
      allocate(turb % v_mean(-nb:nc));  turb % v_mean = 0.
      allocate(turb % w_mean(-nb:nc));  turb % w_mean = 0.
      allocate(turb % p_mean(-nb:nc));  turb % p_mean = 0.
      if(Flow % heat_transfer) then
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
      if(Flow % heat_transfer) then
        allocate(turb % t2_res(-nb:nc));  turb % t2_res = 0.
        allocate(turb % ut_res(-nb:nc));  turb % ut_res = 0.
        allocate(turb % vt_res(-nb:nc));  turb % vt_res = 0.
        allocate(turb % wt_res(-nb:nc));  turb % wt_res = 0.
      end if ! Flow % heat_transfer

    end if ! turb % statistics

  end if ! HYBRID_LES_PRANDTL

  !----------------!
  !   Wale model   !
  !----------------!
  if(turb % model .eq. LES_WALE) then

    allocate(turb % wale_v(-nb:nc));  turb % wale_v = 0.

    ! Other variables such as time scale, length scale and production
    allocate(turb % t_scale (-nb:nc));  turb % t_scale = 0.
    allocate(turb % l_scale (-nb:nc));  turb % l_scale = 0.
    allocate(turb % vis_t   (-nb:nc));  turb % vis_t   = 0.
    allocate(turb % vis_w   (-nb:nc));  turb % vis_w   = 0.
    allocate(turb % p_kin   (-nb:nc));  turb % p_kin   = 0.
    allocate(turb % y_plus  (-nb:nc));  turb % y_plus  = 0.

    if(Flow % heat_transfer) then
      call Var_Mod_Allocate_New_Only(turb % ut, Grid, 'UT')
      call Var_Mod_Allocate_New_Only(turb % vt, Grid, 'VT')
      call Var_Mod_Allocate_New_Only(turb % wt, Grid, 'WT')
      allocate(turb % con_w(-nb:nc));  turb % con_w = 0.  ! wall cond
    end if ! Flow % heat_transfer

    if(turb % statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(turb % u_mean(-nb:nc));  turb % u_mean = 0.
      allocate(turb % v_mean(-nb:nc));  turb % v_mean = 0.
      allocate(turb % w_mean(-nb:nc));  turb % w_mean = 0.
      allocate(turb % p_mean(-nb:nc));  turb % p_mean = 0.
      if(Flow % heat_transfer) then
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
      if(Flow % heat_transfer) then
        allocate(turb % t2_res(-nb:nc));  turb % t2_res = 0.
        allocate(turb % ut_res(-nb:nc));  turb % ut_res = 0.
        allocate(turb % vt_res(-nb:nc));  turb % vt_res = 0.
        allocate(turb % wt_res(-nb:nc));  turb % wt_res = 0.
      end if ! Flow % heat_transfer

    end if ! turb % statistics

  end if ! LES_WALE

  !-------------------!
  !   Dynamic model   !
  !-------------------!
  if(turb % model .eq. LES_DYNAMIC) then

    allocate(turb % c_dyn(-nb:nc));  turb % c_dyn = 0.

    ! Other variables such as time scale, length scale and production
    allocate(turb % t_scale (-nb:nc));  turb % t_scale = 0.
    allocate(turb % l_scale (-nb:nc));  turb % l_scale = 0.
    allocate(turb % vis_t   (-nb:nc));  turb % vis_t   = 0.
    allocate(turb % vis_w   (-nb:nc));  turb % vis_w   = 0.  ! wall visc
    allocate(turb % p_kin   (-nb:nc));  turb % p_kin   = 0.
    allocate(turb % y_plus  (-nb:nc));  turb % y_plus  = 0.

    if(Flow % heat_transfer) then

      call Var_Mod_Allocate_New_Only(turb % ut, Grid, 'UT')
      call Var_Mod_Allocate_New_Only(turb % vt, Grid, 'VT')
      call Var_Mod_Allocate_New_Only(turb % wt, Grid, 'WT')
      allocate(turb % con_w(-nb:nc));  turb % con_w = 0.  ! wall cond

    end if ! Flow % heat_transfer

    if(turb % statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(turb % u_mean(-nb:nc));  turb % u_mean = 0.
      allocate(turb % v_mean(-nb:nc));  turb % v_mean = 0.
      allocate(turb % w_mean(-nb:nc));  turb % w_mean = 0.
      allocate(turb % p_mean(-nb:nc));  turb % p_mean = 0.
      if(Flow % heat_transfer) then
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
      if(Flow % heat_transfer) then
        allocate(turb % t2_res(-nb:nc));  turb % t2_res = 0.
        allocate(turb % ut_res(-nb:nc));  turb % ut_res = 0.
        allocate(turb % vt_res(-nb:nc));  turb % vt_res = 0.
        allocate(turb % wt_res(-nb:nc));  turb % wt_res = 0.
      end if ! Flow % heat_transfer

    end if ! turb % statistics

  end if ! LES_DYNAMIC

  !---------------------------------!
  !   Direct numerical simulation   !
  !---------------------------------!
  if(turb % model .eq. DNS .or.  &
     turb % model .eq. NO_TURBULENCE_MODEL) then

    if(turb % statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(turb % u_mean(-nb:nc));  turb % u_mean = 0.
      allocate(turb % v_mean(-nb:nc));  turb % v_mean = 0.
      allocate(turb % w_mean(-nb:nc));  turb % w_mean = 0.
      allocate(turb % p_mean(-nb:nc));  turb % p_mean = 0.
      if(Flow % heat_transfer) then
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
      if(Flow % heat_transfer) then
        allocate(turb % t2_res(-nb:nc));  turb % t2_res = 0.
        allocate(turb % ut_res(-nb:nc));  turb % ut_res = 0.
        allocate(turb % vt_res(-nb:nc));  turb % vt_res = 0.
        allocate(turb % wt_res(-nb:nc));  turb % wt_res = 0.
      end if ! Flow % heat_transfer

    end if ! turb % statistics

  end if ! DNS

  !-------------------!
  !   Hybrid model    !
  !-------------------!
  if(turb % model .eq. HYBRID_LES_RANS) then

    ! Main model's variables (for RANS part)
    call Var_Mod_Allocate_Solution(turb % kin,  Grid, 'KIN',  '')
    call Var_Mod_Allocate_Solution(turb % eps,  Grid, 'EPS',  '')
    call Var_Mod_Allocate_Solution(turb % zeta, Grid, 'ZETA', '')
    call Var_Mod_Allocate_Solution(turb % f22,  Grid, 'F22',  '')

    ! Main model's variables (for LES part)
    allocate(turb % c_dyn(-nb:nc));  turb % c_dyn = 0.

    ! Other variables such as time scale, length scale and production
    allocate(turb % vis_t_eff(-nb:nc));  turb % vis_t_eff = 0.
    allocate(turb % vis_t_sgs(-nb:nc));  turb % vis_t_sgs = 0.
    allocate(turb % vis_t    (-nb:nc));  turb % vis_t     = 0.
    allocate(turb % vis_w    (-nb:nc));  turb % vis_w     = 0.  ! wall visc
    allocate(turb % t_scale  (-nb:nc));  turb % t_scale   = 0.
    allocate(turb % l_scale  (-nb:nc));  turb % l_scale   = 0.
    allocate(turb % alpha_l  (-nb:nc));  turb % alpha_l   = 0.
    allocate(turb % alpha_u  (-nb:nc));  turb % alpha_u   = 0.
    allocate(turb % p_kin    (-nb:nc));  turb % p_kin     = 0.
    allocate(turb % y_plus   (-nb:nc));  turb % y_plus    = 0.

    ! Hydraulic roughness given by formula
    if(turb % rough_walls) then
      allocate(turb % z_o_f  (-nb:nc)); turb % z_o_f   = -1.0
      allocate(turb % id_zone(-nb:nc)); turb % id_zone =  0.0
    end if

    if(Flow % heat_transfer) then
      call Var_Mod_Allocate_Solution(turb % t2, Grid, 'T2', '')
      call Var_Mod_Allocate_New_Only(turb % ut, Grid, 'UT')
      call Var_Mod_Allocate_New_Only(turb % vt, Grid, 'VT')
      call Var_Mod_Allocate_New_Only(turb % wt, Grid, 'WT')
      allocate(turb % con_w(-nb:nc));  turb % con_w = 0.  ! wall cond
      allocate(turb % p_t2 (-nb:nc));  turb % p_t2  = 0.
    end if
 
    if(Flow % n_scalars > 0.or.Flow % heat_transfer) then            
      ! Reynolds stresses
      call Var_Mod_Allocate_New_Only(turb % uu, Grid, 'UU')
      call Var_Mod_Allocate_New_Only(turb % vv, Grid, 'VV')
      call Var_Mod_Allocate_New_Only(turb % ww, Grid, 'WW')
      call Var_Mod_Allocate_New_Only(turb % uv, Grid, 'UV')
      call Var_Mod_Allocate_New_Only(turb % uw, Grid, 'UW')
      call Var_Mod_Allocate_New_Only(turb % vw, Grid, 'VW')
    end if ! Flow % heat_transfer

    !  Wall difussivity for user scalar
    if(Flow % n_scalars > 0) then            
      allocate(turb % diff_w(-nb:nc));  turb % diff_w = 0.  
      allocate(turb % uc(-nb:nc));  turb % uc = 0.  
      allocate(turb % vc(-nb:nc));  turb % vc = 0.  
      allocate(turb % wc(-nb:nc));  turb % wc = 0.  
    end if


    if(turb % statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(turb % u_mean(-nb:nc));  turb % u_mean = 0.
      allocate(turb % v_mean(-nb:nc));  turb % v_mean = 0.
      allocate(turb % w_mean(-nb:nc));  turb % w_mean = 0.
      allocate(turb % p_mean(-nb:nc));  turb % p_mean = 0.
      if(Flow % heat_transfer) then
        allocate(turb % t_mean(-nb:nc));  turb % t_mean = 0.
      end if

      ! Time-averaged modeled quantities
      allocate(turb % kin_mean (-nb:nc));  turb % kin_mean  = 0.
      allocate(turb % eps_mean (-nb:nc));  turb % eps_mean  = 0.
      allocate(turb % f22_mean (-nb:nc));  turb % f22_mean  = 0.
      allocate(turb % zeta_mean(-nb:nc));  turb % zeta_mean = 0.
      if(Flow % heat_transfer) then
        allocate(turb % t2_mean(-nb:nc));  turb % t2_mean = 0.
      end if

      ! Resolved Reynolds stresses
      allocate(turb % uu_res(-nb:nc));  turb % uu_res = 0.
      allocate(turb % vv_res(-nb:nc));  turb % vv_res = 0.
      allocate(turb % ww_res(-nb:nc));  turb % ww_res = 0.
      allocate(turb % uv_res(-nb:nc));  turb % uv_res = 0.
      allocate(turb % vw_res(-nb:nc));  turb % vw_res = 0.
      allocate(turb % uw_res(-nb:nc));  turb % uw_res = 0.

      ! Resolved turbulent heat fluxes
      if(Flow % heat_transfer) then
        allocate(turb % t2_res(-nb:nc));  turb % t2_res = 0.
        allocate(turb % ut_res(-nb:nc));  turb % ut_res = 0.
        allocate(turb % vt_res(-nb:nc));  turb % vt_res = 0.
        allocate(turb % wt_res(-nb:nc));  turb % wt_res = 0.
        allocate(turb % ut_mean(-nb:nc));  turb % ut_mean = 0.
        allocate(turb % vt_mean(-nb:nc));  turb % vt_mean = 0.
        allocate(turb % wt_mean(-nb:nc));  turb % wt_mean = 0.
      end if ! Flow % heat_transfer

    end if ! turb % statistics

    if(Flow % buoyancy .eq. THERMALLY_DRIVEN) then
      allocate(turb % g_buoy(-nb:nc));  turb % g_buoy = 0.
    end if ! Flow % buoyancy .eq. THERMALLY_DRIVEN

  end if ! HYBRID_LES_RANS

  !----------------------------------------------------------!
  !   Scalars are independent of the turbulence model used   !
  !----------------------------------------------------------!
  if(turb % statistics) then
    if(Flow % n_scalars > 0) then
      allocate(turb % scalar_mean(Flow % n_scalars, -nb:nc))
      turb % scalar_mean = 0.
    end if
  end if

  end subroutine

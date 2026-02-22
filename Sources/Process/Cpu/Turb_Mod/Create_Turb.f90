!==============================================================================!
  subroutine Create_Turb(Turb, Flow)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
  type(Field_Type), target :: Flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Matrix_Type), pointer :: A
  integer                    :: nb, nc
!==============================================================================!

  ! Give some sign
  if(First_Proc()) then
    print '(a)', ' # Creating the turbulence module'
  end if

  ! Store pointers
  Turb % pnt_flow => Flow
  Turb % pnt_grid => Flow % pnt_grid

  ! Take aliases
  Grid => Flow % pnt_grid
  A    => Flow % pnt_matrix
  nb = Grid % n_bnd_cells
  nc = Grid % n_cells

  ! Create deltas
  if(Turb % model .eq. SPALART_ALLMARAS .or.  &
     Turb % model .eq. DES_SPALART      .or.  &
     Turb % model .eq. HYBRID_LES_PRANDTL) then
    allocate(Turb % h_max(-nb:nc));  Turb % h_max = 0.
    allocate(Turb % h_min(-nb:nc));  Turb % h_min = 0.
    allocate(Turb % h_w  (-nb:nc));  Turb % h_w   = 0.
  end if

  allocate(Turb % tau_wall(-nb:nc));  Turb % tau_wall = 0.
  allocate(Turb % y_plus  (-nb:nc));  Turb % y_plus   = 0.

  ! Hydraulic roughness
  allocate(Turb % z_o(-nb:nc)); Turb % z_o = 0.0

  !  Wall difussivity for user scalar
  if(Flow % n_scalars > 0) then
    allocate(Turb % diff_w(-nb:nc));  Turb % diff_w = 0.
    allocate(Turb % uc    (-nb:nc));  Turb % uc = 0.
    allocate(Turb % vc    (-nb:nc));  Turb % vc = 0.
    allocate(Turb % wc    (-nb:nc));  Turb % wc = 0.
  end if

  !-----------------!
  !   K-eps model   !
  !-----------------!
  if(Turb % model .eq. K_EPS) then

    ! Variables we solve for: k and epsilon
    call Var_Mod_Create_Solution(Turb % kin, A, 'KIN', '')
    call Var_Mod_Create_Solution(Turb % eps, A, 'EPS', '')

    ! Other turbulent quantities
    allocate(Turb % vis_t  (-nb:nc));  Turb % vis_t   = 0.
    allocate(Turb % vis_w  (-nb:nc));  Turb % vis_w   = 0.  ! wall visc
    allocate(Turb % p_kin  (-nb:nc));  Turb % p_kin   = 0.
    allocate(Turb % t_scale(-nb:nc));  Turb % t_scale = 0.

    ! Reynolds stresses
    call Var_Mod_Create_New_Only(Turb % uu, Grid, 'UU')
    call Var_Mod_Create_New_Only(Turb % vv, Grid, 'VV')
    call Var_Mod_Create_New_Only(Turb % ww, Grid, 'WW')
    call Var_Mod_Create_New_Only(Turb % uv, Grid, 'UV')
    call Var_Mod_Create_New_Only(Turb % uw, Grid, 'UW')
    call Var_Mod_Create_New_Only(Turb % vw, Grid, 'VW')

    if(Flow % heat_transfer) then
      call Var_Mod_Create_Solution(Turb % t2, A, 'T2', '')
      call Var_Mod_Create_New_Only(Turb % ut, Grid, 'UT')
      call Var_Mod_Create_New_Only(Turb % vt, Grid, 'VT')
      call Var_Mod_Create_New_Only(Turb % wt, Grid, 'WT')
      allocate(Turb % con_w(-nb:nc));  Turb % con_w = 0.  ! wall cond
      allocate(Turb % p_t2 (-nb:nc));  Turb % p_t2  = 0.
    end if ! Flow % heat_transfer

    ! Turbulent statistics; if needed
    if(Turb % statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(Turb % u_mean(-nb:nc));  Turb % u_mean = 0.
      allocate(Turb % v_mean(-nb:nc));  Turb % v_mean = 0.
      allocate(Turb % w_mean(-nb:nc));  Turb % w_mean = 0.
      allocate(Turb % p_mean(-nb:nc));  Turb % p_mean = 0.
      if(Flow % heat_transfer) then
        allocate(Turb % t_mean (-nb:nc));  Turb % t_mean  = 0.
        allocate(Turb % q_mean (-nb:nc));  Turb % q_mean  = 0.
        allocate(Turb % t2_mean(-nb:nc));  Turb % t2_mean = 0.
        allocate(Turb % ut_mean(-nb:nc));  Turb % ut_mean = 0.
        allocate(Turb % vt_mean(-nb:nc));  Turb % vt_mean = 0.
        allocate(Turb % wt_mean(-nb:nc));  Turb % wt_mean = 0.
      end if

      ! Time-averaged modeled quantities
      allocate(Turb % kin_mean(-nb:nc));  Turb % kin_mean = 0.
      allocate(Turb % eps_mean(-nb:nc));  Turb % eps_mean = 0.

      ! Resolved Reynolds stresses
      allocate(Turb % uu_res(-nb:nc));  Turb % uu_res = 0.
      allocate(Turb % vv_res(-nb:nc));  Turb % vv_res = 0.
      allocate(Turb % ww_res(-nb:nc));  Turb % ww_res = 0.
      allocate(Turb % uv_res(-nb:nc));  Turb % uv_res = 0.
      allocate(Turb % vw_res(-nb:nc));  Turb % vw_res = 0.
      allocate(Turb % uw_res(-nb:nc));  Turb % uw_res = 0.

      ! Resolved turbulent heat fluxes
      if(Flow % heat_transfer) then
        allocate(Turb % t2_res (-nb:nc));  Turb % t2_res  = 0.
        allocate(Turb % ut_res (-nb:nc));  Turb % ut_res  = 0.
        allocate(Turb % vt_res (-nb:nc));  Turb % vt_res  = 0.
        allocate(Turb % wt_res (-nb:nc));  Turb % wt_res  = 0.
      end if ! Flow % heat_transfer

    end if ! Turb % statistics

    if(Flow % buoyancy .eq. THERMALLY_DRIVEN) then
      allocate(Turb % g_buoy(-nb:nc));  Turb % g_buoy = 0.
    end if ! Flow % buoyancy .eq. THERMALLY_DRIVEN

  end if ! K_EPS

  !-----------------!
  !   K-omega-sst   !
  !-----------------!
  if(Turb % model .eq. K_OMEGA_SST) then
    call Var_Mod_Create_Solution(Turb % kin,   A, 'KIN', '')
    call Var_Mod_Create_Solution(Turb % omega, A, 'OMG', '')

    ! Other turbulent quantities
    allocate(Turb % vis_t  (-nb:nc));  Turb % vis_t   = 0.
    allocate(Turb % vis_w  (-nb:nc));  Turb % vis_w   = 0.  ! wall visc
    allocate(Turb % p_kin  (-nb:nc));  Turb % p_kin   = 0.
    allocate(Turb % t_scale(-nb:nc));  Turb % t_scale = 0.

    ! Reynolds stresses
    call Var_Mod_Create_New_Only(Turb % uu, Grid, 'UU')
    call Var_Mod_Create_New_Only(Turb % vv, Grid, 'VV')
    call Var_Mod_Create_New_Only(Turb % ww, Grid, 'WW')
    call Var_Mod_Create_New_Only(Turb % uv, Grid, 'UV')
    call Var_Mod_Create_New_Only(Turb % uw, Grid, 'UW')
    call Var_Mod_Create_New_Only(Turb % vw, Grid, 'VW')

    if(Flow % heat_transfer) then
      call Var_Mod_Create_Solution(Turb % t2, A, 'T2', '')
      call Var_Mod_Create_New_Only(Turb % ut, Grid, 'UT')
      call Var_Mod_Create_New_Only(Turb % vt, Grid, 'VT')
      call Var_Mod_Create_New_Only(Turb % wt, Grid, 'WT')
      allocate(Turb % con_w(-nb:nc));  Turb % con_w = 0.
      allocate(Turb % p_t2 (-nb:nc));  Turb % p_t2  = 0.
    end if ! Flow % heat_transfer

    ! Turbulent statistics; if needed
    if(Turb % statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(Turb % u_mean(-nb:nc));  Turb % u_mean = 0.
      allocate(Turb % v_mean(-nb:nc));  Turb % v_mean = 0.
      allocate(Turb % w_mean(-nb:nc));  Turb % w_mean = 0.
      allocate(Turb % p_mean(-nb:nc));  Turb % p_mean = 0.
      if(Flow % heat_transfer) then
        allocate(Turb % t_mean (-nb:nc));  Turb % t_mean  = 0.
        allocate(Turb % q_mean (-nb:nc));  Turb % q_mean  = 0.
        allocate(Turb % t2_mean(-nb:nc));  Turb % t2_mean = 0.
        allocate(Turb % ut_mean(-nb:nc));  Turb % ut_mean = 0.
        allocate(Turb % vt_mean(-nb:nc));  Turb % vt_mean = 0.
        allocate(Turb % wt_mean(-nb:nc));  Turb % wt_mean = 0.
      end if

      ! Time-averaged modeled quantities
      allocate(Turb % kin_mean(-nb:nc))  ;  Turb % kin_mean   = 0.
      allocate(Turb % omega_mean(-nb:nc));  Turb % omega_mean = 0.
    end if
  end if


  !------------------!
  !   K-eps-zeta-f   !
  !------------------!
  if(Turb % model .eq. K_EPS_ZETA_F) then

    ! Main model's variables
    call Var_Mod_Create_Solution(Turb % kin,  A, 'KIN',  '')
    call Var_Mod_Create_Solution(Turb % eps,  A, 'EPS',  '')
    call Var_Mod_Create_Solution(Turb % zeta, A, 'ZETA', '')
    call Var_Mod_Create_Solution(Turb % f22,  A, 'F22',  '')

    ! Other variables such as time scale, length scale and production
    allocate(Turb % t_scale (-nb:nc));  Turb % t_scale = 0.
    allocate(Turb % l_scale (-nb:nc));  Turb % l_scale = 0.
    allocate(Turb % alpha_l (-nb:nc));  Turb % alpha_l = 0.
    allocate(Turb % alpha_u (-nb:nc));  Turb % alpha_u = 0.
    allocate(Turb % vis_t   (-nb:nc));  Turb % vis_t   = 0.
    allocate(Turb % vis_w   (-nb:nc));  Turb % vis_w   = 0.  ! wall visc
    allocate(Turb % p_kin   (-nb:nc));  Turb % p_kin   = 0.

    ! Reynolds stresses
    call Var_Mod_Create_New_Only(Turb % uu, Grid, 'UU')
    call Var_Mod_Create_New_Only(Turb % vv, Grid, 'VV')
    call Var_Mod_Create_New_Only(Turb % ww, Grid, 'WW')
    call Var_Mod_Create_New_Only(Turb % uv, Grid, 'UV')
    call Var_Mod_Create_New_Only(Turb % uw, Grid, 'UW')
    call Var_Mod_Create_New_Only(Turb % vw, Grid, 'VW')

    if(Flow % heat_transfer) then
      call Var_Mod_Create_Solution(Turb % t2, A, 'T2', '')
      call Var_Mod_Create_New_Only(Turb % ut, Grid, 'UT')
      call Var_Mod_Create_New_Only(Turb % vt, Grid, 'VT')
      call Var_Mod_Create_New_Only(Turb % wt, Grid, 'WT')
      allocate(Turb % con_w(-nb:nc));  Turb % con_w = 0.  ! wall cond
      allocate(Turb % p_t2 (-nb:nc));  Turb % p_t2  = 0.
    end if ! Flow % heat_transfer

    if(Turb % statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(Turb % u_mean(-nb:nc));  Turb % u_mean = 0.
      allocate(Turb % v_mean(-nb:nc));  Turb % v_mean = 0.
      allocate(Turb % w_mean(-nb:nc));  Turb % w_mean = 0.
      allocate(Turb % p_mean(-nb:nc));  Turb % p_mean = 0.
      if(Flow % heat_transfer) then
        allocate(Turb % t_mean(-nb:nc));  Turb % t_mean = 0.
        allocate(Turb % q_mean(-nb:nc));  Turb % q_mean = 0.
      end if

      ! Time-averaged modeled quantities
      allocate(Turb % kin_mean (-nb:nc));  Turb % kin_mean  = 0.
      allocate(Turb % eps_mean (-nb:nc));  Turb % eps_mean  = 0.
      allocate(Turb % f22_mean (-nb:nc));  Turb % f22_mean  = 0.
      allocate(Turb % zeta_mean(-nb:nc));  Turb % zeta_mean = 0.
      if(Flow % heat_transfer) then
        allocate(Turb % t2_mean(-nb:nc));  Turb % t2_mean   = 0.
      end if

      ! Resolved Reynolds stresses
      allocate(Turb % uu_res(-nb:nc));  Turb % uu_res = 0.
      allocate(Turb % vv_res(-nb:nc));  Turb % vv_res = 0.
      allocate(Turb % ww_res(-nb:nc));  Turb % ww_res = 0.
      allocate(Turb % uv_res(-nb:nc));  Turb % uv_res = 0.
      allocate(Turb % vw_res(-nb:nc));  Turb % vw_res = 0.
      allocate(Turb % uw_res(-nb:nc));  Turb % uw_res = 0.

      ! Resolved turbulent heat fluxes
      if(Flow % heat_transfer) then
        allocate(Turb % t2_res (-nb:nc));  Turb % t2_res  = 0.
        allocate(Turb % ut_res (-nb:nc));  Turb % ut_res  = 0.
        allocate(Turb % vt_res (-nb:nc));  Turb % vt_res  = 0.
        allocate(Turb % wt_res (-nb:nc));  Turb % wt_res  = 0.
        allocate(Turb % ut_mean(-nb:nc));  Turb % ut_mean = 0.
        allocate(Turb % vt_mean(-nb:nc));  Turb % vt_mean = 0.
        allocate(Turb % wt_mean(-nb:nc));  Turb % wt_mean = 0.
      end if ! Flow % heat_transfer

    end if ! Turb % statistics

    if(Flow % buoyancy .eq. THERMALLY_DRIVEN) then
      allocate(Turb % g_buoy(-nb:nc));  Turb % g_buoy = 0.
    end if ! Flow % buoyancy .eq. THERMALLY_DRIVEN

  end if ! K_EPS_ZETA_F

  !----------------------!
  !   Spalart Allmaras   !
  !----------------------!
  if(Turb % model .eq. SPALART_ALLMARAS .or.  &
     Turb % model .eq. DES_SPALART) then

    call Var_Mod_Create_Solution(Turb % vis, A, 'VIS', '')

    ! Other variables such as time scale, length scale and production
    allocate(Turb % t_scale (-nb:nc));  Turb % t_scale = 0.
    allocate(Turb % l_scale (-nb:nc));  Turb % l_scale = 0.
    allocate(Turb % vis_t   (-nb:nc));  Turb % vis_t   = 0.
    allocate(Turb % vis_w   (-nb:nc));  Turb % vis_w   = 0.  ! wall visc
    allocate(Turb % p_kin   (-nb:nc));  Turb % p_kin   = 0.

    if(Flow % heat_transfer) then
      call Var_Mod_Create_New_Only(Turb % ut, Grid, 'UT')
      call Var_Mod_Create_New_Only(Turb % vt, Grid, 'VT')
      call Var_Mod_Create_New_Only(Turb % wt, Grid, 'WT')
      allocate(Turb % con_w(-nb:nc));  Turb % con_w = 0.  ! wall cond
    end if ! Flow % heat_transfer

    ! Turbulence statistics, if needed
    if(Turb % statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(Turb % u_mean(-nb:nc));  Turb % u_mean = 0.
      allocate(Turb % v_mean(-nb:nc));  Turb % v_mean = 0.
      allocate(Turb % w_mean(-nb:nc));  Turb % w_mean = 0.
      allocate(Turb % p_mean(-nb:nc));  Turb % p_mean = 0.
      if(Flow % heat_transfer) then
        allocate(Turb % t_mean(-nb:nc));  Turb % t_mean = 0.
        allocate(Turb % q_mean(-nb:nc));  Turb % q_mean = 0.
      end if

      ! Time-averaged modeled quantities
      allocate(Turb % vis_mean (-nb:nc));  Turb % vis_mean = 0.

      ! Resolved Reynolds stresses
      allocate(Turb % uu_res(-nb:nc));  Turb % uu_res = 0.
      allocate(Turb % vv_res(-nb:nc));  Turb % vv_res = 0.
      allocate(Turb % ww_res(-nb:nc));  Turb % ww_res = 0.
      allocate(Turb % uv_res(-nb:nc));  Turb % uv_res = 0.
      allocate(Turb % vw_res(-nb:nc));  Turb % vw_res = 0.
      allocate(Turb % uw_res(-nb:nc));  Turb % uw_res = 0.

      ! Resolved turbulent heat fluxes
      if(Flow % heat_transfer) then
        allocate(Turb % t2_res(-nb:nc));  Turb % t2_res = 0.
        allocate(Turb % ut_res(-nb:nc));  Turb % ut_res = 0.
        allocate(Turb % vt_res(-nb:nc));  Turb % vt_res = 0.
        allocate(Turb % wt_res(-nb:nc));  Turb % wt_res = 0.
      end if ! Flow % heat_transfer

    end if ! Turb % statistics

  end if ! SPALART_ALLMARAS & DES_SPALART

  !-----------------------!
  !   Smagorinsky model   !
  !-----------------------!
  if(Turb % model .eq. LES_SMAGORINSKY) then

    ! Variables such as time scale, length scale and production
    allocate(Turb % t_scale (-nb:nc));  Turb % t_scale = 0.
    allocate(Turb % l_scale (-nb:nc));  Turb % l_scale = 0.
    allocate(Turb % vis_t   (-nb:nc));  Turb % vis_t   = 0.
    allocate(Turb % vis_w   (-nb:nc));  Turb % vis_w   = 0.  ! wall visc
    allocate(Turb % p_kin   (-nb:nc));  Turb % p_kin   = 0.

    if(Flow % heat_transfer) then
      call Var_Mod_Create_New_Only(Turb % ut, Grid, 'UT')
      call Var_Mod_Create_New_Only(Turb % vt, Grid, 'VT')
      call Var_Mod_Create_New_Only(Turb % wt, Grid, 'WT')
      allocate(Turb % con_w(-nb:nc));  Turb % con_w = 0.  ! wall cond
    end if ! Flow % heat_transfer

    if(Turb % statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(Turb % u_mean(-nb:nc));  Turb % u_mean = 0.
      allocate(Turb % v_mean(-nb:nc));  Turb % v_mean = 0.
      allocate(Turb % w_mean(-nb:nc));  Turb % w_mean = 0.
      allocate(Turb % p_mean(-nb:nc));  Turb % p_mean = 0.
      if(Flow % heat_transfer) then
        allocate(Turb % t_mean(-nb:nc));  Turb % t_mean = 0.
        allocate(Turb % q_mean(-nb:nc));  Turb % q_mean = 0.
      end if

      ! Resolved Reynolds stresses
      allocate(Turb % uu_res(-nb:nc));  Turb % uu_res = 0.
      allocate(Turb % vv_res(-nb:nc));  Turb % vv_res = 0.
      allocate(Turb % ww_res(-nb:nc));  Turb % ww_res = 0.
      allocate(Turb % uv_res(-nb:nc));  Turb % uv_res = 0.
      allocate(Turb % vw_res(-nb:nc));  Turb % vw_res = 0.
      allocate(Turb % uw_res(-nb:nc));  Turb % uw_res = 0.

      ! Resolved turbulent heat fluxes
      if(Flow % heat_transfer) then
        allocate(Turb % t2_res(-nb:nc));  Turb % t2_res = 0.
        allocate(Turb % ut_res(-nb:nc));  Turb % ut_res = 0.
        allocate(Turb % vt_res(-nb:nc));  Turb % vt_res = 0.
        allocate(Turb % wt_res(-nb:nc));  Turb % wt_res = 0.
      end if ! Flow % heat_transfer

    end if ! Turb % statistics

  end if ! LES_SMAGORINSKY


  !------------------------------!
  !   Hybrid_Les_Prandtl model   !
  !------------------------------!
  if(Turb % model .eq. HYBRID_LES_PRANDTL) then

    ! Dynamic Smagorinsky constant for particle SGS models
    allocate(Turb % c_dyn(-nb:nc));  Turb % c_dyn = 0.

    ! Other variables such as time scale, length scale and production
    allocate(Turb % t_scale (-nb:nc));  Turb % t_scale = 0.
    allocate(Turb % l_scale (-nb:nc));  Turb % l_scale = 0.
    allocate(Turb % vis_t   (-nb:nc));  Turb % vis_t   = 0.
    allocate(Turb % vis_w   (-nb:nc));  Turb % vis_w   = 0.  ! wall visc
    allocate(Turb % p_kin   (-nb:nc));  Turb % p_kin   = 0.

    if(Flow % heat_transfer) then
      call Var_Mod_Create_New_Only(Turb % ut, Grid, 'UT')
      call Var_Mod_Create_New_Only(Turb % vt, Grid, 'VT')
      call Var_Mod_Create_New_Only(Turb % wt, Grid, 'WT')
      allocate(Turb % con_w(-nb:nc));  Turb % con_w = 0.  ! wall cond
    end if ! Flow % heat_transfer

    if(Turb % statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(Turb % u_mean(-nb:nc));  Turb % u_mean = 0.
      allocate(Turb % v_mean(-nb:nc));  Turb % v_mean = 0.
      allocate(Turb % w_mean(-nb:nc));  Turb % w_mean = 0.
      allocate(Turb % p_mean(-nb:nc));  Turb % p_mean = 0.
      if(Flow % heat_transfer) then
        allocate(Turb % t_mean(-nb:nc));  Turb % t_mean = 0.
      end if

      ! Resolved Reynolds stresses
      allocate(Turb % uu_res(-nb:nc));  Turb % uu_res = 0.
      allocate(Turb % vv_res(-nb:nc));  Turb % vv_res = 0.
      allocate(Turb % ww_res(-nb:nc));  Turb % ww_res = 0.
      allocate(Turb % uv_res(-nb:nc));  Turb % uv_res = 0.
      allocate(Turb % vw_res(-nb:nc));  Turb % vw_res = 0.
      allocate(Turb % uw_res(-nb:nc));  Turb % uw_res = 0.

      ! Resolved turbulent heat fluxes
      if(Flow % heat_transfer) then
        allocate(Turb % t2_res(-nb:nc));  Turb % t2_res = 0.
        allocate(Turb % ut_res(-nb:nc));  Turb % ut_res = 0.
        allocate(Turb % vt_res(-nb:nc));  Turb % vt_res = 0.
        allocate(Turb % wt_res(-nb:nc));  Turb % wt_res = 0.
      end if ! Flow % heat_transfer

    end if ! Turb % statistics

  end if ! HYBRID_LES_PRANDTL

  !----------------!
  !   Wale model   !
  !----------------!
  if(Turb % model .eq. LES_WALE) then

    allocate(Turb % wale_v(-nb:nc));  Turb % wale_v = 0.

    ! Other variables such as time scale, length scale and production
    allocate(Turb % t_scale (-nb:nc));  Turb % t_scale = 0.
    allocate(Turb % l_scale (-nb:nc));  Turb % l_scale = 0.
    allocate(Turb % vis_t   (-nb:nc));  Turb % vis_t   = 0.
    allocate(Turb % vis_w   (-nb:nc));  Turb % vis_w   = 0.
    allocate(Turb % p_kin   (-nb:nc));  Turb % p_kin   = 0.

    if(Flow % heat_transfer) then
      call Var_Mod_Create_New_Only(Turb % ut, Grid, 'UT')
      call Var_Mod_Create_New_Only(Turb % vt, Grid, 'VT')
      call Var_Mod_Create_New_Only(Turb % wt, Grid, 'WT')
      allocate(Turb % con_w(-nb:nc));  Turb % con_w = 0.  ! wall cond
    end if ! Flow % heat_transfer

    if(Turb % statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(Turb % u_mean(-nb:nc));  Turb % u_mean = 0.
      allocate(Turb % v_mean(-nb:nc));  Turb % v_mean = 0.
      allocate(Turb % w_mean(-nb:nc));  Turb % w_mean = 0.
      allocate(Turb % p_mean(-nb:nc));  Turb % p_mean = 0.
      if(Flow % heat_transfer) then
        allocate(Turb % t_mean(-nb:nc));  Turb % t_mean = 0.
        allocate(Turb % q_mean(-nb:nc));  Turb % q_mean = 0.
      end if

      ! Resolved Reynolds stresses
      allocate(Turb % uu_res(-nb:nc));  Turb % uu_res = 0.
      allocate(Turb % vv_res(-nb:nc));  Turb % vv_res = 0.
      allocate(Turb % ww_res(-nb:nc));  Turb % ww_res = 0.
      allocate(Turb % uv_res(-nb:nc));  Turb % uv_res = 0.
      allocate(Turb % vw_res(-nb:nc));  Turb % vw_res = 0.
      allocate(Turb % uw_res(-nb:nc));  Turb % uw_res = 0.

      ! Resolved turbulent heat fluxes
      if(Flow % heat_transfer) then
        allocate(Turb % t2_res(-nb:nc));  Turb % t2_res = 0.
        allocate(Turb % ut_res(-nb:nc));  Turb % ut_res = 0.
        allocate(Turb % vt_res(-nb:nc));  Turb % vt_res = 0.
        allocate(Turb % wt_res(-nb:nc));  Turb % wt_res = 0.
      end if ! Flow % heat_transfer

    end if ! Turb % statistics

  end if ! LES_WALE

  !-------------------!
  !   Dynamic model   !
  !-------------------!
  if(Turb % model .eq. LES_DYNAMIC) then

    allocate(Turb % c_dyn(-nb:nc));  Turb % c_dyn = 0.

    ! Other variables such as time scale, length scale and production
    allocate(Turb % t_scale (-nb:nc));  Turb % t_scale = 0.
    allocate(Turb % l_scale (-nb:nc));  Turb % l_scale = 0.
    allocate(Turb % vis_t   (-nb:nc));  Turb % vis_t   = 0.
    allocate(Turb % vis_w   (-nb:nc));  Turb % vis_w   = 0.  ! wall visc
    allocate(Turb % p_kin   (-nb:nc));  Turb % p_kin   = 0.

    if(Flow % heat_transfer) then

      call Var_Mod_Create_New_Only(Turb % ut, Grid, 'UT')
      call Var_Mod_Create_New_Only(Turb % vt, Grid, 'VT')
      call Var_Mod_Create_New_Only(Turb % wt, Grid, 'WT')
      allocate(Turb % con_w(-nb:nc));  Turb % con_w = 0.  ! wall cond

    end if ! Flow % heat_transfer

    if(Turb % statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(Turb % u_mean(-nb:nc));  Turb % u_mean = 0.
      allocate(Turb % v_mean(-nb:nc));  Turb % v_mean = 0.
      allocate(Turb % w_mean(-nb:nc));  Turb % w_mean = 0.
      allocate(Turb % p_mean(-nb:nc));  Turb % p_mean = 0.
      if(Flow % heat_transfer) then
        allocate(Turb % t_mean(-nb:nc));  Turb % t_mean = 0.
        allocate(Turb % q_mean(-nb:nc));  Turb % q_mean = 0.
      end if

      ! Resolved Reynolds stresses
      allocate(Turb % uu_res(-nb:nc));  Turb % uu_res = 0.
      allocate(Turb % vv_res(-nb:nc));  Turb % vv_res = 0.
      allocate(Turb % ww_res(-nb:nc));  Turb % ww_res = 0.
      allocate(Turb % uv_res(-nb:nc));  Turb % uv_res = 0.
      allocate(Turb % vw_res(-nb:nc));  Turb % vw_res = 0.
      allocate(Turb % uw_res(-nb:nc));  Turb % uw_res = 0.

      ! Resolved turbulent heat fluxes
      if(Flow % heat_transfer) then
        allocate(Turb % t2_res(-nb:nc));  Turb % t2_res = 0.
        allocate(Turb % ut_res(-nb:nc));  Turb % ut_res = 0.
        allocate(Turb % vt_res(-nb:nc));  Turb % vt_res = 0.
        allocate(Turb % wt_res(-nb:nc));  Turb % wt_res = 0.
      end if ! Flow % heat_transfer

    end if ! Turb % statistics

  end if ! LES_DYNAMIC

  !-------------------------------!
  !   Tensorial viscosity model   !
  !-------------------------------!
  if(Turb % model .eq. LES_TVM) then

    ! Variables such as time scale, length scale and production
    allocate(Turb % t_scale (-nb:nc));  Turb % t_scale = 0.
    allocate(Turb % l_scale (-nb:nc));  Turb % l_scale = 0.
    allocate(Turb % vis_t   (-nb:nc));  Turb % vis_t   = 0.
    allocate(Turb % vis_w   (-nb:nc));  Turb % vis_w   = 0.  ! wall visc
    allocate(Turb % p_kin   (-nb:nc));  Turb % p_kin   = 0.

    ! Tensorial turbulent viscosity
    allocate(Turb % ten_turb_11 (-nb:nc));  Turb % ten_turb_11 = 0.
    allocate(Turb % ten_turb_12 (-nb:nc));  Turb % ten_turb_12 = 0.
    allocate(Turb % ten_turb_13 (-nb:nc));  Turb % ten_turb_13 = 0.
    allocate(Turb % ten_turb_21 (-nb:nc));  Turb % ten_turb_21 = 0.
    allocate(Turb % ten_turb_22 (-nb:nc));  Turb % ten_turb_22 = 0.
    allocate(Turb % ten_turb_23 (-nb:nc));  Turb % ten_turb_23 = 0.
    allocate(Turb % ten_turb_31 (-nb:nc));  Turb % ten_turb_31 = 0.
    allocate(Turb % ten_turb_32 (-nb:nc));  Turb % ten_turb_32 = 0.
    allocate(Turb % ten_turb_33 (-nb:nc));  Turb % ten_turb_33 = 0.

    ! Turbulent stress tensor
    allocate(Turb % tau_11      (-nb:nc));  Turb % tau_11      = 0.
    allocate(Turb % tau_12      (-nb:nc));  Turb % tau_12      = 0.
    allocate(Turb % tau_13      (-nb:nc));  Turb % tau_13      = 0.
    allocate(Turb % tau_21      (-nb:nc));  Turb % tau_21      = 0.
    allocate(Turb % tau_22      (-nb:nc));  Turb % tau_22      = 0.
    allocate(Turb % tau_23      (-nb:nc));  Turb % tau_23      = 0.
    allocate(Turb % tau_31      (-nb:nc));  Turb % tau_31      = 0.
    allocate(Turb % tau_32      (-nb:nc));  Turb % tau_32      = 0.
    allocate(Turb % tau_33      (-nb:nc));  Turb % tau_33      = 0.

    if(Flow % heat_transfer) then

      call Var_Mod_Create_New_Only(Turb % ut, Grid, 'UT')
      call Var_Mod_Create_New_Only(Turb % vt, Grid, 'VT')
      call Var_Mod_Create_New_Only(Turb % wt, Grid, 'WT')
      allocate(Turb % con_w(-nb:nc));  Turb % con_w = 0.  ! wall cond

    end if ! Flow % heat_transfer

    if(Turb % statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(Turb % u_mean(-nb:nc));  Turb % u_mean = 0.
      allocate(Turb % v_mean(-nb:nc));  Turb % v_mean = 0.
      allocate(Turb % w_mean(-nb:nc));  Turb % w_mean = 0.
      allocate(Turb % p_mean(-nb:nc));  Turb % p_mean = 0.
      if(Flow % heat_transfer) then
        allocate(Turb % t_mean(-nb:nc));  Turb % t_mean = 0.
        allocate(Turb % q_mean(-nb:nc));  Turb % q_mean = 0.
      end if

      ! Resolved Reynolds stresses
      allocate(Turb % uu_res(-nb:nc));  Turb % uu_res = 0.
      allocate(Turb % vv_res(-nb:nc));  Turb % vv_res = 0.
      allocate(Turb % ww_res(-nb:nc));  Turb % ww_res = 0.
      allocate(Turb % uv_res(-nb:nc));  Turb % uv_res = 0.
      allocate(Turb % vw_res(-nb:nc));  Turb % vw_res = 0.
      allocate(Turb % uw_res(-nb:nc));  Turb % uw_res = 0.

      ! Resolved turbulent heat fluxes
      if(Flow % heat_transfer) then
        allocate(Turb % t2_res(-nb:nc));  Turb % t2_res = 0.
        allocate(Turb % ut_res(-nb:nc));  Turb % ut_res = 0.
        allocate(Turb % vt_res(-nb:nc));  Turb % vt_res = 0.
        allocate(Turb % wt_res(-nb:nc));  Turb % wt_res = 0.
      end if ! Flow % heat_transfer

    end if ! Turb % statistics

  end if ! LES_TVM

  !---------------------------------!
  !   Direct numerical simulation   !
  !---------------------------------!
  if(Turb % model .eq. DNS .or.  &
     Turb % model .eq. NO_TURBULENCE_MODEL) then

    if(Turb % statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(Turb % u_mean(-nb:nc));  Turb % u_mean = 0.
      allocate(Turb % v_mean(-nb:nc));  Turb % v_mean = 0.
      allocate(Turb % w_mean(-nb:nc));  Turb % w_mean = 0.
      allocate(Turb % p_mean(-nb:nc));  Turb % p_mean = 0.
      if(Flow % heat_transfer) then
        allocate(Turb % t_mean(-nb:nc));  Turb % t_mean = 0.
        allocate(Turb % q_mean(-nb:nc));  Turb % q_mean = 0.
      end if

      ! Resolved Reynolds stresses
      allocate(Turb % uu_res(-nb:nc));  Turb % uu_res = 0.
      allocate(Turb % vv_res(-nb:nc));  Turb % vv_res = 0.
      allocate(Turb % ww_res(-nb:nc));  Turb % ww_res = 0.
      allocate(Turb % uv_res(-nb:nc));  Turb % uv_res = 0.
      allocate(Turb % vw_res(-nb:nc));  Turb % vw_res = 0.
      allocate(Turb % uw_res(-nb:nc));  Turb % uw_res = 0.

      ! Resolved turbulent heat fluxes
      if(Flow % heat_transfer) then
        allocate(Turb % t2_res(-nb:nc));  Turb % t2_res = 0.
        allocate(Turb % ut_res(-nb:nc));  Turb % ut_res = 0.
        allocate(Turb % vt_res(-nb:nc));  Turb % vt_res = 0.
        allocate(Turb % wt_res(-nb:nc));  Turb % wt_res = 0.
      end if ! Flow % heat_transfer

    end if ! Turb % statistics

  end if ! DNS

  !-------------------!
  !   Hybrid model    !
  !-------------------!
  if(Turb % model .eq. HYBRID_LES_RANS) then

    ! Main model's variables (for RANS part)
    call Var_Mod_Create_Solution(Turb % kin,  A, 'KIN',  '')
    call Var_Mod_Create_Solution(Turb % eps,  A, 'EPS',  '')
    call Var_Mod_Create_Solution(Turb % zeta, A, 'ZETA', '')
    call Var_Mod_Create_Solution(Turb % f22,  A, 'F22',  '')

    ! Main model's variables (for LES part)
    allocate(Turb % c_dyn(-nb:nc));  Turb % c_dyn = 0.

    ! Other variables such as time scale, length scale and production
    allocate(Turb % vis_t_eff(-nb:nc));  Turb % vis_t_eff = 0.
    allocate(Turb % vis_t_sgs(-nb:nc));  Turb % vis_t_sgs = 0.
    allocate(Turb % vis_t    (-nb:nc));  Turb % vis_t     = 0.
    allocate(Turb % vis_w    (-nb:nc));  Turb % vis_w     = 0.  ! wall visc
    allocate(Turb % t_scale  (-nb:nc));  Turb % t_scale   = 0.
    allocate(Turb % l_scale  (-nb:nc));  Turb % l_scale   = 0.
    allocate(Turb % alpha_l  (-nb:nc));  Turb % alpha_l   = 0.
    allocate(Turb % alpha_u  (-nb:nc));  Turb % alpha_u   = 0.
    allocate(Turb % p_kin    (-nb:nc));  Turb % p_kin     = 0.

    if(Flow % heat_transfer) then
      call Var_Mod_Create_Solution(Turb % t2, A, 'T2', '')
      call Var_Mod_Create_New_Only(Turb % ut, Grid, 'UT')
      call Var_Mod_Create_New_Only(Turb % vt, Grid, 'VT')
      call Var_Mod_Create_New_Only(Turb % wt, Grid, 'WT')
      allocate(Turb % con_w(-nb:nc));  Turb % con_w = 0.  ! wall cond
      allocate(Turb % p_t2 (-nb:nc));  Turb % p_t2  = 0.
    end if

    if(Flow % n_scalars > 0.or.Flow % heat_transfer) then
      ! Reynolds stresses
      call Var_Mod_Create_New_Only(Turb % uu, Grid, 'UU')
      call Var_Mod_Create_New_Only(Turb % vv, Grid, 'VV')
      call Var_Mod_Create_New_Only(Turb % ww, Grid, 'WW')
      call Var_Mod_Create_New_Only(Turb % uv, Grid, 'UV')
      call Var_Mod_Create_New_Only(Turb % uw, Grid, 'UW')
      call Var_Mod_Create_New_Only(Turb % vw, Grid, 'VW')
    end if ! Flow % heat_transfer



    if(Turb % statistics) then

      ! Time-averaged velocities (and temperature)
      allocate(Turb % u_mean(-nb:nc));  Turb % u_mean = 0.
      allocate(Turb % v_mean(-nb:nc));  Turb % v_mean = 0.
      allocate(Turb % w_mean(-nb:nc));  Turb % w_mean = 0.
      allocate(Turb % p_mean(-nb:nc));  Turb % p_mean = 0.
      if(Flow % heat_transfer) then
        allocate(Turb % t_mean(-nb:nc));  Turb % t_mean = 0.
        allocate(Turb % q_mean(-nb:nc));  Turb % q_mean = 0.
      end if

      ! Time-averaged modeled quantities
      allocate(Turb % kin_mean (-nb:nc));  Turb % kin_mean  = 0.
      allocate(Turb % eps_mean (-nb:nc));  Turb % eps_mean  = 0.
      allocate(Turb % f22_mean (-nb:nc));  Turb % f22_mean  = 0.
      allocate(Turb % zeta_mean(-nb:nc));  Turb % zeta_mean = 0.
      if(Flow % heat_transfer) then
        allocate(Turb % t2_mean(-nb:nc));  Turb % t2_mean = 0.
      end if

      ! Resolved Reynolds stresses
      allocate(Turb % uu_res(-nb:nc));  Turb % uu_res = 0.
      allocate(Turb % vv_res(-nb:nc));  Turb % vv_res = 0.
      allocate(Turb % ww_res(-nb:nc));  Turb % ww_res = 0.
      allocate(Turb % uv_res(-nb:nc));  Turb % uv_res = 0.
      allocate(Turb % vw_res(-nb:nc));  Turb % vw_res = 0.
      allocate(Turb % uw_res(-nb:nc));  Turb % uw_res = 0.

      ! Resolved turbulent heat fluxes
      if(Flow % heat_transfer) then
        allocate(Turb % t2_res (-nb:nc));  Turb % t2_res  = 0.
        allocate(Turb % ut_res (-nb:nc));  Turb % ut_res  = 0.
        allocate(Turb % vt_res (-nb:nc));  Turb % vt_res  = 0.
        allocate(Turb % wt_res (-nb:nc));  Turb % wt_res  = 0.
        allocate(Turb % ut_mean(-nb:nc));  Turb % ut_mean = 0.
        allocate(Turb % vt_mean(-nb:nc));  Turb % vt_mean = 0.
        allocate(Turb % wt_mean(-nb:nc));  Turb % wt_mean = 0.
      end if ! Flow % heat_transfer

    end if ! Turb % statistics

    if(Flow % buoyancy .eq. THERMALLY_DRIVEN) then
      allocate(Turb % g_buoy(-nb:nc));  Turb % g_buoy = 0.
    end if ! Flow % buoyancy .eq. THERMALLY_DRIVEN

  end if ! HYBRID_LES_RANS

  !----------------------------------------------------------!
  !   Scalars are independent of the turbulence model used   !
  !----------------------------------------------------------!
  if(Turb % statistics) then
    if(Flow % n_scalars > 0) then
      allocate(Turb % scalar_mean(Flow % n_scalars, -nb:nc))
      Turb % scalar_mean = 0.
    end if
  end if

  end subroutine

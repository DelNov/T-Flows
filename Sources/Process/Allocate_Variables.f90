!==============================================================================!
  subroutine Allocate_Variables(grid)
!------------------------------------------------------------------------------!
!   Allocates memory for variables. It is called either from LoaRes            !
!   or from Processor.                                                         !
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Les_Mod
  use Rans_Mod
  use Grid_Mod
  use Grad_Mod
  use User_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)   :: grid
!==============================================================================!

  ! Allocate memory for velocity components ...
  call Var_Mod_Allocate_Solution('U', u, grid)
  call Var_Mod_Allocate_Solution('V', v, grid)
  call Var_Mod_Allocate_Solution('W', w, grid)

  ! ... and their gradients
  call Var_Mod_Allocate_Gradients(u)
  call Var_Mod_Allocate_Gradients(v)
  call Var_Mod_Allocate_Gradients(w)

  ! Allocate memory for pressure correction and pressure
  call Var_Mod_Allocate_New_Only('P',  p,  grid)
  call Var_Mod_Allocate_New_Only('PP', pp, grid)

  ! Pressure gradients are needed too
  call Var_Mod_Allocate_Gradients(p)

  ! It is always calling this - probably not needed
  call Var_Mod_Allocate_Statistics(u)
  call Var_Mod_Allocate_Statistics(v)
  call Var_Mod_Allocate_Statistics(w)
  call Var_Mod_Allocate_Statistics(p)

  allocate(phi_face(grid % n_faces)); phi_face = 0.

  allocate(phi_max(-grid % n_bnd_cells:grid % n_cells)); phi_max = 0.
  allocate(phi_min(-grid % n_bnd_cells:grid % n_cells)); phi_min = 0.

  allocate(flux(grid % n_faces));  flux = 0.

  call Grad_Mod_Allocate_Memory(grid)

  allocate(nearest_wall_cell(-grid % n_bnd_cells:grid % n_cells))
  nearest_wall_cell = 0

  allocate(vis_wall(-grid % n_bnd_cells:grid % n_cells)); vis_wall = 0.

  ! For solution of temperature
  call Control_Mod_Heat_Transfer(verbose = .true.)
  if(heat_transfer .eq. YES) then
    call Var_Mod_Allocate_Solution('T',  t,  grid)
    call Var_Mod_Allocate_Solution('UT', ut, grid)
    call Var_Mod_Allocate_Solution('VT', vt, grid)
    call Var_Mod_Allocate_Solution('WT', wt, grid)
    allocate(con_wall(-grid % n_bnd_cells:grid % n_cells)); con_wall = 0.
  end if

  call Control_Mod_Buoyancy(verbose = .true.)

  call Control_Mod_Turbulence_Model(verbose = .true.)
  call Control_Mod_Turbulence_Model_Variant(verbose = .true.)

  ! Allocate wall distance for all models
  if(turbulence_model .ne. NONE) then
    allocate(y_plus(-grid % n_bnd_cells:grid % n_cells));  y_plus = 0.
  end if

  !----------------------------!
  !   Reynolds stress models   !
  !----------------------------!
  if(turbulence_model .eq. REYNOLDS_STRESS .or.  &
     turbulence_model .eq. HANJALIC_JAKIRLIC) then
    if(turbulence_model_variant .eq. URANS) then
!
!     Should something be done here?
!
    end if
    if(turbulence_model .eq. HANJALIC_JAKIRLIC) then
      allocate(eps_tot(-grid % n_bnd_cells:grid % n_cells)); eps_tot = 0.
    end if

    ! Reynolds stresses
    call Var_Mod_Allocate_Solution('UU', uu, grid)
    call Var_Mod_Allocate_Solution('VV', vv, grid)
    call Var_Mod_Allocate_Solution('WW', ww, grid)
    call Var_Mod_Allocate_Solution('UV', uv, grid)
    call Var_Mod_Allocate_Solution('UW', uw, grid)
    call Var_Mod_Allocate_Solution('VW', vw, grid)

    call Var_Mod_Allocate_New_Only('KIN', kin, grid)
    call Var_Mod_Allocate_Solution('EPS', eps, grid)

    ! Time scale, length scale and production
    allocate(t_scale(-grid % n_bnd_cells:grid % n_cells));  t_scale = 0.
    allocate(l_scale(-grid % n_bnd_cells:grid % n_cells));  l_scale = 0.
    allocate(p_kin  (-grid % n_bnd_cells:grid % n_cells));  p_kin   = 0.

    if(turbulence_model .eq. REYNOLDS_STRESS) then
      call Var_Mod_Allocate_Solution('F22', f22, grid)
      call Var_Mod_Allocate_Gradients(f22)
    else
      call Var_Mod_Allocate_New_Only('F22', f22, grid)
    end if

    if(turbulence_model_variant .eq. URANS) then
      call Var_Mod_Allocate_Statistics(uu)
      call Var_Mod_Allocate_Statistics(vv)
      call Var_Mod_Allocate_Statistics(ww)
      call Var_Mod_Allocate_Statistics(uv)
      call Var_Mod_Allocate_Statistics(uw)
      call Var_Mod_Allocate_Statistics(vw)
      call Var_Mod_Allocate_Statistics(kin)
    end if
  end if  ! turbulence_model .eq. 'EBM' or 'HJ'

  ! Variables for Rans models
  if(turbulence_model .eq. K_EPS) then
    call Var_Mod_Allocate_Solution('KIN', kin, grid)
    call Var_Mod_Allocate_Solution('EPS', eps, grid)

    allocate(u_tau     (-grid % n_bnd_cells:grid % n_cells));  u_tau      = 0.
    allocate(u_tau_mean(-grid % n_bnd_cells:grid % n_cells));  u_tau_mean = 0.
    allocate(p_kin     (-grid % n_bnd_cells:grid % n_cells));  p_kin      = 0.

    if(turbulence_model_variant .eq. URANS) then
      allocate(kin % mean(grid % n_cells));   kin % mean = 0.
      allocate(eps % mean(grid % n_cells));   eps % mean = 0.
    end if
  end if

  if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
     turbulence_model .eq. HYBRID_K_EPS_ZETA_F) then
    call Var_Mod_Allocate_Solution('KIN',  kin,  grid)
    call Var_Mod_Allocate_Solution('EPS',  eps,  grid)
    call Var_Mod_Allocate_Solution('ZETA', zeta, grid)
    call Var_Mod_Allocate_Solution('F22',  f22,  grid)

    allocate(t_scale   (-grid % n_bnd_cells:grid % n_cells));  t_scale    = 0.
    allocate(l_scale   (-grid % n_bnd_cells:grid % n_cells));  l_scale    = 0.
    allocate(u_tau     (-grid % n_bnd_cells:grid % n_cells));  u_tau      = 0.
    allocate(u_tau_mean(-grid % n_bnd_cells:grid % n_cells));  u_tau_mean = 0.
    allocate(p_kin     (-grid % n_bnd_cells:grid % n_cells));  p_kin      = 0.

    if(turbulence_model_variant .eq. URANS) then
      allocate(kin  % mean(grid % n_cells));  kin  % mean = 0.
      allocate(eps  % mean(grid % n_cells));  eps  % mean = 0.
      allocate(zeta % mean(grid % n_cells));  zeta % mean = 0.
      allocate(f22  % mean(grid % n_cells));  f22  % mean = 0.
    end if

    if(buoyancy .eq. YES) then
      call Var_Mod_Allocate_Solution('TT', tt, grid)
      call Var_Mod_Allocate_Statistics(tt)
      allocate(g_buoy   (-grid % n_bnd_cells:grid % n_cells));  g_buoy     = 0.
      allocate(buoy_beta(-grid % n_bnd_cells:grid % n_cells));  buoy_beta  = 0.
      allocate(p_buoy   (-grid % n_bnd_cells:grid % n_cells));  p_buoy     = 0.
      allocate(kin%mean (-grid % n_bnd_cells:grid % n_cells));  kin % mean = 0.
      allocate(eps%mean (-grid % n_bnd_cells:grid % n_cells));  eps % mean = 0.
    end if
  end if

  if(turbulence_model .eq. HYBRID_K_EPS_ZETA_F) then
    allocate(vis_t_sgs(-grid % n_bnd_cells:grid % n_cells));  vis_t_sgs = 0.
    allocate(vis_t_eff(-grid % n_bnd_cells:grid % n_cells));  vis_t_eff = 0.
  end if

  if(turbulence_model .eq. DES_SPALART) then
    allocate(kin_sgs(-grid % n_bnd_cells:grid % n_cells));  kin_sgs= 0.
  end if

  if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
     turbulence_model .eq. DES_SPALART) then
    call Var_Mod_Allocate_Solution('VIS', vis, grid)
  end if

  if(turbulence_model .eq. DES_SPALART) then
    allocate(VIS % mean(grid % n_cells));   VIS % mean= 0.
  end if

  ! Variables defined in Les_Mod.h90:
  if(turbulence_model .eq. LES .or.  &
     turbulence_model .eq. HYBRID_K_EPS_ZETA_F) then
    if(turbulence_model_variant .eq. WALE) then
      allocate(wale_v(-grid % n_bnd_cells:grid % n_cells));  wale_v = 0.
    end if
    if(turbulence_model_variant .eq. DYNAMIC) then
      allocate(u % filt(-grid % n_bnd_cells:grid % n_cells));  u % filt = 0.
      allocate(v % filt(-grid % n_bnd_cells:grid % n_cells));  v % filt = 0.
      allocate(w % filt(-grid % n_bnd_cells:grid % n_cells));  w % filt = 0.

      allocate(c_dyn(-grid % n_bnd_cells:grid % n_cells)); c_dyn = 0.
    end if
    allocate(kin_sgs   (-grid % n_bnd_cells:grid % n_cells)); kin_sgs= 0.
    allocate(c_dyn_mean(-grid % n_bnd_cells:grid % n_cells)); c_dyn_mean = 0.
  end if

  if(turbulence_model .eq. LES .or.  &
     turbulence_model .eq. DNS .or.  &
     turbulence_model .eq. DES_SPALART) then
    allocate(uu % mean(-grid % n_bnd_cells:grid % n_cells)); uu % mean= 0.
    allocate(vv % mean(-grid % n_bnd_cells:grid % n_cells)); vv % mean= 0.
    allocate(ww % mean(-grid % n_bnd_cells:grid % n_cells)); ww % mean= 0.
    allocate(uv % mean(-grid % n_bnd_cells:grid % n_cells)); uv % mean= 0.
    allocate(uw % mean(-grid % n_bnd_cells:grid % n_cells)); uw % mean= 0.
    allocate(vw % mean(-grid % n_bnd_cells:grid % n_cells)); vw % mean= 0.

    allocate(vis_t_mean(grid % n_cells)); vis_t_mean = 0.
    allocate(shear_mean(grid % n_cells)); shear_mean = 0.

    if(heat_transfer .eq. YES) then
      allocate(t %  mean(-grid % n_bnd_cells:grid % n_cells)); t  % mean = 0.
      allocate(tt % mean(-grid % n_bnd_cells:grid % n_cells)); tt % mean = 0.
      allocate(ut % mean(-grid % n_bnd_cells:grid % n_cells)); ut % mean = 0.
      allocate(vt % mean(-grid % n_bnd_cells:grid % n_cells)); vt % mean = 0.
      allocate(wt % mean(-grid % n_bnd_cells:grid % n_cells)); wt % mean = 0.
    end if
  end if

  if(turbulence_model .eq. HYBRID_K_EPS_ZETA_F) then
    allocate(uu % mean(-grid % n_bnd_cells:grid % n_cells)); uu % mean = 0.
    allocate(vv % mean(-grid % n_bnd_cells:grid % n_cells)); vv % mean = 0.
    allocate(ww % mean(-grid % n_bnd_cells:grid % n_cells)); ww % mean = 0.
    allocate(uv % mean(-grid % n_bnd_cells:grid % n_cells)); uv % mean = 0.
    allocate(uw % mean(-grid % n_bnd_cells:grid % n_cells)); uw % mean = 0.
    allocate(vw % mean(-grid % n_bnd_cells:grid % n_cells)); vw % mean = 0.

    allocate(vis_t_mean(grid % n_cells));  vis_t_mean = 0.
    allocate(shear_mean(grid % n_cells));  shear_mean = 0.

    if(heat_transfer .eq. YES) then
      allocate(t  % mean(-grid % n_bnd_cells:grid % n_cells)); t  % mean = 0.
      allocate(tt % mean(-grid % n_bnd_cells:grid % n_cells)); tt % mean = 0.
      allocate(ut % mean(-grid % n_bnd_cells:grid % n_cells)); ut % mean = 0.
      allocate(vt % mean(-grid % n_bnd_cells:grid % n_cells)); vt % mean = 0.
      allocate(wt % mean(-grid % n_bnd_cells:grid % n_cells)); wt % mean = 0.
    end if
  end if

  allocate(vis_t   (-grid % n_bnd_cells:grid % n_cells));  vis_t    = 0.
  allocate(vort    (-grid % n_bnd_cells:grid % n_cells));  vort     = 0.
  allocate(shear   (-grid % n_bnd_cells:grid % n_cells));  shear    = 0.
  allocate(tau_wall(grid % n_cells));                      tau_wall = 0.

  ! User scalars and arrays
  call Control_Mod_Number_Of_User_Scalars(n_user_scalars, verbose = .true.)
  call Control_Mod_Number_Of_User_Arrays (n_user_arrays,  verbose = .true.)
  call User_Mod_Allocate(grid)

  end subroutine

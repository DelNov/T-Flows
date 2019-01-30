!==============================================================================!
  subroutine Turbulence_Allocate(flow)
!------------------------------------------------------------------------------!
!   Allocates memory for variables. It is called either from LoaRes            !
!   or from Processor.                                                         !
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Field_Mod, only: Field_Type, heat_transfer, buoyancy
  use Les_Mod
  use Rans_Mod
  use Grid_Mod,  only: Grid_Type
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, v, w, t, p
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  u    => flow % u
  v    => flow % v
  w    => flow % w
  t    => flow % t
  p    => flow % p

  call Control_Mod_Buoyancy        (verbose = .true.)
  call Control_Mod_Turbulence_Model(verbose = .true.)

  !-----------------!
  !   K-eps model   !
  !-----------------!
  if(turbulence_model .eq. K_EPS) then

    ! Variables we solve for: k and epsilon
    call Var_Mod_Allocate_Solution('KIN', kin, grid)
    call Var_Mod_Allocate_Solution('EPS', eps, grid)

    ! Other turbulent quantities
    allocate(u_tau   (-grid % n_bnd_cells:grid % n_cells));  u_tau    = 0.
    allocate(p_kin   (-grid % n_bnd_cells:grid % n_cells));  p_kin    = 0.
    allocate(vort    (-grid % n_bnd_cells:grid % n_cells));  vort     = 0.
    allocate(shear   (-grid % n_bnd_cells:grid % n_cells));  shear    = 0.
    allocate(vis_t   (-grid % n_bnd_cells:grid % n_cells));  vis_t    = 0.
    allocate(vis_wall(-grid % n_bnd_cells:grid % n_cells));  vis_wall = 0.
    allocate(tau_wall(grid % n_cells));                      tau_wall = 0.
    allocate(y_plus  (-grid % n_bnd_cells:grid % n_cells));  y_plus   = 0.

    if(heat_transfer) then

      call Var_Mod_Allocate_New_Only('UT', ut, grid)
      call Var_Mod_Allocate_New_Only('VT', vt, grid)
      call Var_Mod_Allocate_New_Only('WT', wt, grid)
      allocate(con_wall(-grid % n_bnd_cells:grid % n_cells)); con_wall = 0.

    end if ! heat_transfer

    ! Turbulent statistics; if needed
    if(turbulence_statistics) then

      ! Variables that we solve for
      call Var_Mod_Allocate_Statistics(kin)
      call Var_Mod_Allocate_Statistics(eps)

      ! First moments
      call Var_Mod_Allocate_Statistics(u)
      call Var_Mod_Allocate_Statistics(v)
      call Var_Mod_Allocate_Statistics(w)
      call Var_Mod_Allocate_Statistics(p)

      ! Second moments
      call Var_Mod_Allocate_New_Only('UU', uu, grid)
      call Var_Mod_Allocate_New_Only('VV', vv, grid)
      call Var_Mod_Allocate_New_Only('WW', ww, grid)
      call Var_Mod_Allocate_New_Only('UV', uv, grid)
      call Var_Mod_Allocate_New_Only('UW', uw, grid)
      call Var_Mod_Allocate_New_Only('VW', vw, grid)
      call Var_Mod_Allocate_Statistics(uu)
      call Var_Mod_Allocate_Statistics(vv)
      call Var_Mod_Allocate_Statistics(ww)
      call Var_Mod_Allocate_Statistics(uv)
      call Var_Mod_Allocate_Statistics(uw)
      call Var_Mod_Allocate_Statistics(vw)

      ! Heat transfer & statistics
      if(heat_transfer) then

        call Var_Mod_Allocate_Statistics(t)   ! new value allocated above
        call Var_Mod_Allocate_Statistics(ut)  ! new value allocated above
        call Var_Mod_Allocate_Statistics(vt)  ! new value allocated above
        call Var_Mod_Allocate_Statistics(wt)  ! new value allocated above
        call Var_Mod_Allocate_New_Only('TT', tt,  grid)
        call Var_Mod_Allocate_Statistics(tt)

      end if ! heat_transfer

    end if ! turbulence_statistics

  end if ! K_EPS

  !------------------!
  !   K-eps-zeta-f   !
  !------------------!
  if(turbulence_model .eq. K_EPS_ZETA_F) then

    ! Main model's variables
    call Var_Mod_Allocate_Solution('KIN',  kin,  grid)
    call Var_Mod_Allocate_Solution('EPS',  eps,  grid)
    call Var_Mod_Allocate_Solution('ZETA', zeta, grid)
    call Var_Mod_Allocate_Solution('F22',  f22,  grid)

    ! Other variables such as time scale, length scale and production
    allocate(t_scale (-grid % n_bnd_cells:grid % n_cells));  t_scale  = 0.
    allocate(l_scale (-grid % n_bnd_cells:grid % n_cells));  l_scale  = 0.
    allocate(u_tau   (-grid % n_bnd_cells:grid % n_cells));  u_tau    = 0.
    allocate(p_kin   (-grid % n_bnd_cells:grid % n_cells));  p_kin    = 0.
    allocate(vort    (-grid % n_bnd_cells:grid % n_cells));  vort     = 0.
    allocate(shear   (-grid % n_bnd_cells:grid % n_cells));  shear    = 0.
    allocate(vis_t   (-grid % n_bnd_cells:grid % n_cells));  vis_t    = 0.
    allocate(vis_wall(-grid % n_bnd_cells:grid % n_cells));  vis_wall = 0.
    allocate(tau_wall(grid % n_cells));                      tau_wall = 0.
    allocate(y_plus  (-grid % n_bnd_cells:grid % n_cells));  y_plus   = 0.

    if(heat_transfer) then

      call Var_Mod_Allocate_Solution('T2', t2, grid)
      call Var_Mod_Allocate_New_Only('UT', ut, grid)
      call Var_Mod_Allocate_New_Only('VT', vt, grid)
      call Var_Mod_Allocate_New_Only('WT', wt, grid)
      allocate(con_wall(-grid % n_bnd_cells:grid % n_cells)); con_wall = 0.
      allocate(p_t2    (-grid % n_bnd_cells:grid % n_cells)); p_t2     = 0.

    end if ! heat_transfer

    if(turbulence_statistics) then

      ! First moments
      call Var_Mod_Allocate_Statistics(u)
      call Var_Mod_Allocate_Statistics(v)
      call Var_Mod_Allocate_Statistics(w)
      call Var_Mod_Allocate_Statistics(p)
      call Var_Mod_Allocate_Statistics(kin)
      call Var_Mod_Allocate_Statistics(eps)
      call Var_Mod_Allocate_Statistics(zeta)
      call Var_Mod_Allocate_Statistics(f22)

      ! Second moments
      call Var_Mod_Allocate_New_Only('UU', uu, grid)
      call Var_Mod_Allocate_New_Only('VV', vv, grid)
      call Var_Mod_Allocate_New_Only('WW', ww, grid)
      call Var_Mod_Allocate_New_Only('UV', uv, grid)
      call Var_Mod_Allocate_New_Only('UW', uw, grid)
      call Var_Mod_Allocate_New_Only('VW', vw, grid)
      call Var_Mod_Allocate_Statistics(uu)
      call Var_Mod_Allocate_Statistics(vv)
      call Var_Mod_Allocate_Statistics(ww)
      call Var_Mod_Allocate_Statistics(uv)
      call Var_Mod_Allocate_Statistics(uw)
      call Var_Mod_Allocate_Statistics(vw)

      ! Heat transfer & statistics
      if(heat_transfer) then

        call Var_Mod_Allocate_Statistics(t)   ! new value allocated above
        call Var_Mod_Allocate_Statistics(ut)  ! new value allocated above
        call Var_Mod_Allocate_Statistics(vt)  ! new value allocated above
        call Var_Mod_Allocate_Statistics(wt)  ! new value allocated above
        call Var_Mod_Allocate_New_Only('TT', tt,  grid)
        call Var_Mod_Allocate_Statistics(tt)
        call Var_Mod_Allocate_Statistics(t2)

      end if ! heat_transfer

    end if ! turbulence_statistics

    if(buoyancy) then

     allocate(g_buoy   (-grid % n_bnd_cells:grid % n_cells));  g_buoy    = 0.
     allocate(buoy_beta(-grid % n_bnd_cells:grid % n_cells));  buoy_beta = 0.
     allocate(g_kin    (-grid % n_bnd_cells:grid % n_cells));  g_kin     = 0.

    end if ! buoyancy

  end if ! K_EPS_ZETA_F

  !----------------------------!
  !   Reynolds stress models   !
  !----------------------------!
  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then

    ! Other variables such as time scale, length scale and production
    allocate(t_scale (-grid % n_bnd_cells:grid % n_cells));  t_scale  = 0.
    allocate(l_scale (-grid % n_bnd_cells:grid % n_cells));  l_scale  = 0.
    allocate(p_kin   (-grid % n_bnd_cells:grid % n_cells));  p_kin    = 0.
    allocate(u_tau   (-grid % n_bnd_cells:grid % n_cells));  u_tau    = 0.
    allocate(vort    (-grid % n_bnd_cells:grid % n_cells));  vort     = 0.
    allocate(shear   (-grid % n_bnd_cells:grid % n_cells));  shear    = 0.
    allocate(vis_t   (-grid % n_bnd_cells:grid % n_cells));  vis_t    = 0.
    allocate(vis_wall(-grid % n_bnd_cells:grid % n_cells));  vis_wall = 0.
    allocate(tau_wall(grid % n_cells));                      tau_wall = 0.
    allocate(y_plus  (-grid % n_bnd_cells:grid % n_cells));  y_plus   = 0.

    ! Reynolds stresses
    call Var_Mod_Allocate_Solution('UU', uu, grid)
    call Var_Mod_Allocate_Solution('VV', vv, grid)
    call Var_Mod_Allocate_Solution('WW', ww, grid)
    call Var_Mod_Allocate_Solution('UV', uv, grid)
    call Var_Mod_Allocate_Solution('UW', uw, grid)
    call Var_Mod_Allocate_Solution('VW', vw, grid)

    call Var_Mod_Allocate_New_Only('KIN', kin, grid)
    call Var_Mod_Allocate_Solution('EPS', eps, grid)

    if(heat_transfer) then

      call Var_Mod_Allocate_New_Only('UT', ut, grid)
      call Var_Mod_Allocate_New_Only('VT', vt, grid)
      call Var_Mod_Allocate_New_Only('WT', wt, grid)
      allocate(con_wall(-grid % n_bnd_cells:grid % n_cells)); con_wall = 0.

    end if ! heat_transfer

    if(turbulence_statistics) then

      call Var_Mod_Allocate_Statistics(uu)
      call Var_Mod_Allocate_Statistics(vv)
      call Var_Mod_Allocate_Statistics(ww)
      call Var_Mod_Allocate_Statistics(uv)
      call Var_Mod_Allocate_Statistics(uw)
      call Var_Mod_Allocate_Statistics(vw)
      call Var_Mod_Allocate_Statistics(kin)
      call Var_Mod_Allocate_Statistics(eps)

      if(heat_transfer) then

        call Var_Mod_Allocate_Statistics(t)   ! new value allocated above
        call Var_Mod_Allocate_Statistics(ut)  ! new value allocated above
        call Var_Mod_Allocate_Statistics(vt)  ! new value allocated above
        call Var_Mod_Allocate_Statistics(wt)  ! new value allocated above
        call Var_Mod_Allocate_New_Only('TT', tt,  grid)
        call Var_Mod_Allocate_Statistics(tt)

      end if ! heat_transfer

    end if ! turbulence_statistics

    if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then

      call Var_Mod_Allocate_Solution('F22', f22, grid)

      if(turbulence_statistics) then
        call Var_Mod_Allocate_Statistics(f22)
      end if ! turbulence_statistics

    end if ! RSM_MANCEAU_HANJALIC

    if(turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
      allocate(eps_tot(-grid % n_bnd_cells:grid % n_cells)); eps_tot = 0.
    end if ! RSM_HANJALIC_JAKIRLIC

  end if ! RSM_MANCEAU_HANJALIC & RSM_HANJALIC_JAKIRLIC

  !----------------------!
  !   Spalart Allmaras   !
  !----------------------!
  if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
     turbulence_model .eq. DES_SPALART) then

    call Var_Mod_Allocate_Solution('VIS', vis, grid)

    ! Other variables such as time scale, length scale and production
    allocate(t_scale (-grid % n_bnd_cells:grid % n_cells));  t_scale  = 0.
    allocate(l_scale (-grid % n_bnd_cells:grid % n_cells));  l_scale  = 0.
    allocate(u_tau   (-grid % n_bnd_cells:grid % n_cells));  u_tau    = 0.
    allocate(p_kin   (-grid % n_bnd_cells:grid % n_cells));  p_kin    = 0.
    allocate(vort    (-grid % n_bnd_cells:grid % n_cells));  vort     = 0.
    allocate(shear   (-grid % n_bnd_cells:grid % n_cells));  shear    = 0.
    allocate(vis_t   (-grid % n_bnd_cells:grid % n_cells));  vis_t    = 0.
    allocate(vis_wall(-grid % n_bnd_cells:grid % n_cells));  vis_wall = 0.
    allocate(tau_wall(grid % n_cells));                      tau_wall = 0.
    allocate(y_plus  (-grid % n_bnd_cells:grid % n_cells));  y_plus   = 0.

    if(heat_transfer) then

      call Var_Mod_Allocate_New_Only('UT', ut, grid)
      call Var_Mod_Allocate_New_Only('VT', vt, grid)
      call Var_Mod_Allocate_New_Only('WT', wt, grid)
      allocate(con_wall(-grid % n_bnd_cells:grid % n_cells)); con_wall = 0.

    end if ! heat_transfer

    ! Turbulence statistics, if needed
    if(turbulence_statistics) then

      call Var_Mod_Allocate_Statistics(vis)

      ! First moments
      call Var_Mod_Allocate_Statistics(u)
      call Var_Mod_Allocate_Statistics(v)
      call Var_Mod_Allocate_Statistics(w)
      call Var_Mod_Allocate_Statistics(p)

      ! Second moments
      call Var_Mod_Allocate_New_Only('UU', uu, grid)
      call Var_Mod_Allocate_New_Only('VV', vv, grid)
      call Var_Mod_Allocate_New_Only('WW', ww, grid)
      call Var_Mod_Allocate_New_Only('UV', uv, grid)
      call Var_Mod_Allocate_New_Only('UW', uw, grid)
      call Var_Mod_Allocate_New_Only('VW', vw, grid)
      call Var_Mod_Allocate_Statistics(uu)
      call Var_Mod_Allocate_Statistics(vv)
      call Var_Mod_Allocate_Statistics(ww)
      call Var_Mod_Allocate_Statistics(uv)
      call Var_Mod_Allocate_Statistics(uw)
      call Var_Mod_Allocate_Statistics(vw)

      ! Heat transfer & statistics
      if(heat_transfer) then

        call Var_Mod_Allocate_Statistics(t)   ! new value allocated above
        call Var_Mod_Allocate_Statistics(ut)  ! new value allocated above
        call Var_Mod_Allocate_Statistics(vt)  ! new value allocated above
        call Var_Mod_Allocate_Statistics(wt)  ! new value allocated above
        call Var_Mod_Allocate_New_Only('TT', tt,  grid)
        call Var_Mod_Allocate_Statistics(tt)

      end if ! heat_transfer

    end if ! turbulence_statistics

  end if ! SPALART_ALLMARAS & DES_SPALART

  !-----------------------!
  !   Smagorinsky model   !
  !-----------------------!
  if(turbulence_model .eq. LES_SMAGORINSKY) then

    allocate(nearest_wall_cell(-grid % n_bnd_cells:grid % n_cells))
    nearest_wall_cell = 0

    ! Other variables such as time scale, length scale and production
    allocate(t_scale (-grid % n_bnd_cells:grid % n_cells));  t_scale  = 0.
    allocate(l_scale (-grid % n_bnd_cells:grid % n_cells));  l_scale  = 0.
    allocate(u_tau   (-grid % n_bnd_cells:grid % n_cells));  u_tau    = 0.
    allocate(p_kin   (-grid % n_bnd_cells:grid % n_cells));  p_kin    = 0.
    allocate(vort    (-grid % n_bnd_cells:grid % n_cells));  vort     = 0.
    allocate(shear   (-grid % n_bnd_cells:grid % n_cells));  shear    = 0.
    allocate(vis_t   (-grid % n_bnd_cells:grid % n_cells));  vis_t    = 0.
    allocate(vis_wall(-grid % n_bnd_cells:grid % n_cells));  vis_wall = 0.
    allocate(tau_wall(grid % n_cells));                      tau_wall = 0.
    allocate(y_plus  (-grid % n_bnd_cells:grid % n_cells));  y_plus   = 0.

    if(heat_transfer) then

      call Var_Mod_Allocate_New_Only('UT', ut, grid)
      call Var_Mod_Allocate_New_Only('VT', vt, grid)
      call Var_Mod_Allocate_New_Only('WT', wt, grid)
      allocate(con_wall(-grid % n_bnd_cells:grid % n_cells)); con_wall = 0.

    end if ! heat_transfer

    if(turbulence_statistics) then

      ! First moments
      call Var_Mod_Allocate_Statistics(u)
      call Var_Mod_Allocate_Statistics(v)
      call Var_Mod_Allocate_Statistics(w)
      call Var_Mod_Allocate_Statistics(p)
      ! Second moments
      call Var_Mod_Allocate_New_Only('UU', uu, grid)
      call Var_Mod_Allocate_New_Only('VV', vv, grid)
      call Var_Mod_Allocate_New_Only('WW', ww, grid)
      call Var_Mod_Allocate_New_Only('UV', uv, grid)
      call Var_Mod_Allocate_New_Only('UW', uw, grid)
      call Var_Mod_Allocate_New_Only('VW', vw, grid)
      call Var_Mod_Allocate_Statistics(uu)
      call Var_Mod_Allocate_Statistics(vv)
      call Var_Mod_Allocate_Statistics(ww)
      call Var_Mod_Allocate_Statistics(uv)
      call Var_Mod_Allocate_Statistics(uw)
      call Var_Mod_Allocate_Statistics(vw)

      ! Heat transfer & statistics
      if(heat_transfer) then

        call Var_Mod_Allocate_Statistics(t)   ! new value allocated above
        call Var_Mod_Allocate_Statistics(ut)  ! new value allocated above
        call Var_Mod_Allocate_Statistics(vt)  ! new value allocated above
        call Var_Mod_Allocate_Statistics(wt)  ! new value allocated above
        call Var_Mod_Allocate_New_Only('TT', tt,  grid)
        call Var_Mod_Allocate_Statistics(tt)

      end if ! heat_transfer

    end if ! turbulence_statistics

  end if ! LES_SMAGORINSKY

  !----------------!
  !   Wale model   !
  !----------------!
  if(turbulence_model .eq. LES_WALE) then

    allocate(wale_v(-grid % n_bnd_cells:grid % n_cells));  wale_v = 0.

    ! Other variables such as time scale, length scale and production
    allocate(t_scale (-grid % n_bnd_cells:grid % n_cells));  t_scale  = 0.
    allocate(l_scale (-grid % n_bnd_cells:grid % n_cells));  l_scale  = 0.
    allocate(u_tau   (-grid % n_bnd_cells:grid % n_cells));  u_tau    = 0.
    allocate(p_kin   (-grid % n_bnd_cells:grid % n_cells));  p_kin    = 0.
    allocate(vort    (-grid % n_bnd_cells:grid % n_cells));  vort     = 0.
    allocate(shear   (-grid % n_bnd_cells:grid % n_cells));  shear    = 0.
    allocate(vis_t   (-grid % n_bnd_cells:grid % n_cells));  vis_t    = 0.
    allocate(vis_wall(-grid % n_bnd_cells:grid % n_cells));  vis_wall = 0.
    allocate(tau_wall(grid % n_cells));                      tau_wall = 0.
    allocate(y_plus  (-grid % n_bnd_cells:grid % n_cells));  y_plus   = 0.

    if(heat_transfer) then

      call Var_Mod_Allocate_New_Only('UT', ut, grid)
      call Var_Mod_Allocate_New_Only('VT', vt, grid)
      call Var_Mod_Allocate_New_Only('WT', wt, grid)
      allocate(con_wall(-grid % n_bnd_cells:grid % n_cells)); con_wall = 0.

    end if ! heat_transfer

    if(turbulence_statistics) then

      ! First moments
      call Var_Mod_Allocate_Statistics(u)
      call Var_Mod_Allocate_Statistics(v)
      call Var_Mod_Allocate_Statistics(w)
      call Var_Mod_Allocate_Statistics(p)

      ! Second moments
      call Var_Mod_Allocate_New_Only('UU', uu, grid)
      call Var_Mod_Allocate_New_Only('VV', vv, grid)
      call Var_Mod_Allocate_New_Only('WW', ww, grid)
      call Var_Mod_Allocate_New_Only('UV', uv, grid)
      call Var_Mod_Allocate_New_Only('UW', uw, grid)
      call Var_Mod_Allocate_New_Only('VW', vw, grid)
      call Var_Mod_Allocate_Statistics(uu)
      call Var_Mod_Allocate_Statistics(vv)
      call Var_Mod_Allocate_Statistics(ww)
      call Var_Mod_Allocate_Statistics(uv)
      call Var_Mod_Allocate_Statistics(uw)
      call Var_Mod_Allocate_Statistics(vw)

      ! Heat transfer & statistics
      if(heat_transfer) then

        call Var_Mod_Allocate_Statistics(t)   ! new value allocated above
        call Var_Mod_Allocate_Statistics(ut)  ! new value allocated above
        call Var_Mod_Allocate_Statistics(vt)  ! new value allocated above
        call Var_Mod_Allocate_Statistics(wt)  ! new value allocated above
        call Var_Mod_Allocate_New_Only('TT', tt,  grid)
        call Var_Mod_Allocate_Statistics(tt)

      end if ! heat_transfer

    end if ! turbulence_statistics

  end if ! LES_WALE

  !-------------------!
  !   Dynamic model   !
  !-------------------!
  if(turbulence_model .eq. LES_DYNAMIC) then

    allocate(c_dyn(-grid % n_bnd_cells:grid % n_cells));  c_dyn = 0.

    ! Other variables such as time scale, length scale and production
    allocate(t_scale (-grid % n_bnd_cells:grid % n_cells));  t_scale  = 0.
    allocate(l_scale (-grid % n_bnd_cells:grid % n_cells));  l_scale  = 0.
    allocate(u_tau   (-grid % n_bnd_cells:grid % n_cells));  u_tau    = 0.
    allocate(p_kin   (-grid % n_bnd_cells:grid % n_cells));  p_kin    = 0.
    allocate(vort    (-grid % n_bnd_cells:grid % n_cells));  vort     = 0.
    allocate(shear   (-grid % n_bnd_cells:grid % n_cells));  shear    = 0.
    allocate(vis_t   (-grid % n_bnd_cells:grid % n_cells));  vis_t    = 0.
    allocate(vis_wall(-grid % n_bnd_cells:grid % n_cells));  vis_wall = 0.
    allocate(tau_wall(grid % n_cells));                      tau_wall = 0.
    allocate(y_plus  (-grid % n_bnd_cells:grid % n_cells));  y_plus   = 0.

    if(heat_transfer) then

      call Var_Mod_Allocate_New_Only('UT', ut, grid)
      call Var_Mod_Allocate_New_Only('VT', vt, grid)
      call Var_Mod_Allocate_New_Only('WT', wt, grid)
      allocate(con_wall(-grid % n_bnd_cells:grid % n_cells)); con_wall = 0.

    end if ! heat_transfer

    if(turbulence_statistics) then

      ! First moments
      call Var_Mod_Allocate_Statistics(u)
      call Var_Mod_Allocate_Statistics(v)
      call Var_Mod_Allocate_Statistics(w)
      call Var_Mod_Allocate_Statistics(p)

      ! Second moments
      call Var_Mod_Allocate_New_Only('UU', uu, grid)
      call Var_Mod_Allocate_New_Only('VV', vv, grid)
      call Var_Mod_Allocate_New_Only('WW', ww, grid)
      call Var_Mod_Allocate_New_Only('UV', uv, grid)
      call Var_Mod_Allocate_New_Only('UW', uw, grid)
      call Var_Mod_Allocate_New_Only('VW', vw, grid)
      call Var_Mod_Allocate_Statistics(uu)
      call Var_Mod_Allocate_Statistics(vv)
      call Var_Mod_Allocate_Statistics(ww)
      call Var_Mod_Allocate_Statistics(uv)
      call Var_Mod_Allocate_Statistics(uw)
      call Var_Mod_Allocate_Statistics(vw)

      ! Heat transfer & statistics
      if(heat_transfer) then

        call Var_Mod_Allocate_Statistics(t)   ! new value allocated above
        call Var_Mod_Allocate_Statistics(ut)  ! new value allocated above
        call Var_Mod_Allocate_Statistics(vt)  ! new value allocated above
        call Var_Mod_Allocate_Statistics(wt)  ! new value allocated above
        call Var_Mod_Allocate_New_Only('TT', tt,  grid)
        call Var_Mod_Allocate_Statistics(tt)

      end if ! heat_transfer

    end if ! turbulence_statistics

  end if ! LES_DYNAMIC

  !-------------------!
  !   Hybrid model    !
  !-------------------!
  if(turbulence_model .eq. HYBRID_LES_RANS) then

    ! Main model's variables (for RANS part)
    call Var_Mod_Allocate_Solution('KIN',  kin,  grid)
    call Var_Mod_Allocate_Solution('EPS',  eps,  grid)
    call Var_Mod_Allocate_Solution('ZETA', zeta, grid)
    call Var_Mod_Allocate_Solution('F22',  f22,  grid)
    call Var_Mod_Allocate_Statistics(kin)
    call Var_Mod_Allocate_Statistics(eps)
    call Var_Mod_Allocate_Statistics(zeta)
    call Var_Mod_Allocate_Statistics(f22)

    ! Main model's variables (for LES part)
    allocate(c_dyn(-grid % n_bnd_cells:grid % n_cells));  c_dyn = 0.

    ! Other variables such as time scale, length scale and production
    allocate(vis_t_eff(-grid % n_bnd_cells:grid % n_cells));  vis_t_eff = 0.
    allocate(vis_t_sgs(-grid % n_bnd_cells:grid % n_cells));  vis_t_sgs = 0.
    allocate(t_scale  (-grid % n_bnd_cells:grid % n_cells));  t_scale   = 0.
    allocate(l_scale  (-grid % n_bnd_cells:grid % n_cells));  l_scale   = 0.
    allocate(u_tau    (-grid % n_bnd_cells:grid % n_cells));  u_tau     = 0.
    allocate(p_kin    (-grid % n_bnd_cells:grid % n_cells));  p_kin     = 0.
    allocate(vort     (-grid % n_bnd_cells:grid % n_cells));  vort      = 0.
    allocate(shear    (-grid % n_bnd_cells:grid % n_cells));  shear     = 0.
    allocate(vis_t    (-grid % n_bnd_cells:grid % n_cells));  vis_t     = 0.
    allocate(vis_wall (-grid % n_bnd_cells:grid % n_cells));  vis_wall  = 0.
    allocate(tau_wall (grid % n_cells));                      tau_wall  = 0.
    allocate(y_plus   (-grid % n_bnd_cells:grid % n_cells));  y_plus    = 0.

    if(heat_transfer) then

      call Var_Mod_Allocate_New_Only('UT', ut, grid)
      call Var_Mod_Allocate_New_Only('VT', vt, grid)
      call Var_Mod_Allocate_New_Only('WT', wt, grid)
      allocate(con_wall(-grid % n_bnd_cells:grid % n_cells)); con_wall = 0.

    end if ! heat_transfer

    if(turbulence_statistics) then

      ! First moments
      call Var_Mod_Allocate_Statistics(u)
      call Var_Mod_Allocate_Statistics(v)
      call Var_Mod_Allocate_Statistics(w)
      call Var_Mod_Allocate_Statistics(p)

      ! Second moments
      call Var_Mod_Allocate_New_Only('UU', uu, grid)
      call Var_Mod_Allocate_New_Only('VV', vv, grid)
      call Var_Mod_Allocate_New_Only('WW', ww, grid)
      call Var_Mod_Allocate_New_Only('UV', uv, grid)
      call Var_Mod_Allocate_New_Only('UW', uw, grid)
      call Var_Mod_Allocate_New_Only('VW', vw, grid)
      call Var_Mod_Allocate_Statistics(uu)
      call Var_Mod_Allocate_Statistics(vv)
      call Var_Mod_Allocate_Statistics(ww)
      call Var_Mod_Allocate_Statistics(uv)
      call Var_Mod_Allocate_Statistics(uw)
      call Var_Mod_Allocate_Statistics(vw)

      ! Heat transfer & statistics
      if(heat_transfer) then

        call Var_Mod_Allocate_Statistics(t)   ! new value allocated above
        call Var_Mod_Allocate_Statistics(ut)  ! new value allocated above
        call Var_Mod_Allocate_Statistics(vt)  ! new value allocated above
        call Var_Mod_Allocate_Statistics(wt)  ! new value allocated above
        call Var_Mod_Allocate_New_Only('TT', tt,  grid)
        call Var_Mod_Allocate_Statistics(tt)

      end if ! heat_transfer

    end if ! turbulence_statistics

    if(buoyancy) then

      call Var_Mod_Allocate_Solution('TT', tt, grid)
      call Var_Mod_Allocate_Statistics(tt)
      allocate(g_buoy   (-grid % n_bnd_cells:grid % n_cells));  g_buoy    = 0.
      allocate(buoy_beta(-grid % n_bnd_cells:grid % n_cells));  buoy_beta = 0.
      allocate(g_kin    (-grid % n_bnd_cells:grid % n_cells));  g_kin     = 0.

    end if ! buoyancy

  end if ! HYBRID_LES_RANS

  end subroutine

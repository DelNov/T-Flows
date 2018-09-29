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

  call Control_Mod_Buoyancy                 (verbose = .true.)
  call Control_Mod_Turbulence_Statistics    (verbose = .true.)
  call Control_Mod_Turbulence_Model         (verbose = .true.)
  call Control_Mod_Turbulence_Wall_Treatment(verbose = .true.)

  ! Gradient matrices are always needed
  call Grad_Mod_Allocate_Memory(grid)

  !----------------------------!
  !                            !
  !   Navier-Stokes equation   !
  !                            !
  !----------------------------!

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

  ! Mass flow rates at cell faces are always needed
  allocate(flux(grid % n_faces));  flux = 0.

  !-----------------------------------------!
  !                                         !
  !   Enthalpy conservation (temperature)   !
  !                                         !
  !-----------------------------------------!
  call Control_Mod_Heat_Transfer(verbose = .true.)
  if(heat_transfer) then
    call Var_Mod_Allocate_Solution('T',  t,  grid)
    allocate(con_wall(-grid % n_bnd_cells:grid % n_cells)); con_wall = 0.

    call Var_Mod_Allocate_New_Only('TT', tt, grid)
    call Var_Mod_Allocate_New_Only('UT', ut, grid)
    call Var_Mod_Allocate_New_Only('VT', vt, grid)
    call Var_Mod_Allocate_New_Only('WT', wt, grid)

    if(turbulence_statistics) then
      call Var_Mod_Allocate_Statistics(t)
      call Var_Mod_Allocate_Statistics(tt)
      call Var_Mod_Allocate_Statistics(ut)
      call Var_Mod_Allocate_Statistics(vt)
      call Var_Mod_Allocate_Statistics(wt)
    end if
  end if


  !-----------------------!
  !                       !
  !   Turbulence models   !
  !                       !
  !-----------------------!

  ! Arrayes needed by all models, or almost all
  if(turbulence_model .ne. NONE) then
    allocate(vort    (-grid % n_bnd_cells:grid % n_cells));  vort     = 0.
    allocate(shear   (-grid % n_bnd_cells:grid % n_cells));  shear    = 0.
    allocate(vis_t   (-grid % n_bnd_cells:grid % n_cells));  vis_t    = 0.
    allocate(vis_wall(-grid % n_bnd_cells:grid % n_cells));  vis_wall = 0.
    allocate(tau_wall(grid % n_cells));                      tau_wall = 0.
    allocate(y_plus  (-grid % n_bnd_cells:grid % n_cells));  y_plus   = 0.
  end if

  !-----------------!
  !   K-eps model   !
  !-----------------!
  if(turbulence_model .eq. K_EPS) then

    ! K and epsilon
    call Var_Mod_Allocate_Solution('KIN', kin, grid)
    call Var_Mod_Allocate_Solution('EPS', eps, grid)

    ! Other turbulent quantities
    allocate(u_tau(-grid % n_bnd_cells:grid % n_cells));  u_tau      = 0.
    allocate(p_kin(-grid % n_bnd_cells:grid % n_cells));  p_kin      = 0.

    ! Turbulent statistics; if needed
    if(turbulence_statistics) then
      call Var_Mod_Allocate_Statistics(kin)               
      call Var_Mod_Allocate_Statistics(eps)               
    end if
  end if

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
    allocate(t_scale(-grid % n_bnd_cells:grid % n_cells));  t_scale = 0.
    allocate(l_scale(-grid % n_bnd_cells:grid % n_cells));  l_scale = 0.
    allocate(u_tau  (-grid % n_bnd_cells:grid % n_cells));  u_tau   = 0.
    allocate(p_kin  (-grid % n_bnd_cells:grid % n_cells));  p_kin   = 0.

    ! Turbulent statistics; if needed
    if(turbulence_statistics) then
      call Var_Mod_Allocate_Statistics(kin)
      call Var_Mod_Allocate_Statistics(eps)
      call Var_Mod_Allocate_Statistics(zeta)
      call Var_Mod_Allocate_Statistics(f22)
      allocate(vis_t_eff(-grid % n_bnd_cells:grid % n_cells));  vis_t_eff = 0.
      allocate(vis_t_sgs(-grid % n_bnd_cells:grid % n_cells));  vis_t_sgs = 0.
    end if

    if(buoyancy) then
      call Var_Mod_Allocate_Solution('TT', tt, grid)
      call Var_Mod_Allocate_Statistics(tt)
      allocate(g_buoy   (-grid % n_bnd_cells:grid % n_cells));  g_buoy    = 0.
      allocate(buoy_beta(-grid % n_bnd_cells:grid % n_cells));  buoy_beta = 0.
      allocate(p_buoy   (-grid % n_bnd_cells:grid % n_cells));  p_buoy    = 0.
    end if

  end if

  !----------------------------!
  !   Reynolds stress models   !
  !----------------------------!
  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then

    ! Reynolds stresses
    call Var_Mod_Allocate_Solution('UU', uu, grid)
    call Var_Mod_Allocate_Solution('VV', vv, grid)
    call Var_Mod_Allocate_Solution('WW', ww, grid)
    call Var_Mod_Allocate_Solution('UV', uv, grid)
    call Var_Mod_Allocate_Solution('UW', uw, grid)
    call Var_Mod_Allocate_Solution('VW', vw, grid)

    call Var_Mod_Allocate_New_Only('KIN', kin, grid)
    call Var_Mod_Allocate_Solution('EPS', eps, grid)

    if(turbulence_statistics) then
      call Var_Mod_Allocate_Statistics(uu)
      call Var_Mod_Allocate_Statistics(vv)
      call Var_Mod_Allocate_Statistics(ww)
      call Var_Mod_Allocate_Statistics(uv)
      call Var_Mod_Allocate_Statistics(uw)
      call Var_Mod_Allocate_Statistics(vw)
      call Var_Mod_Allocate_Statistics(kin)
      call Var_Mod_Allocate_Statistics(eps)
    end if

    if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
      call Var_Mod_Allocate_Solution('F22', f22, grid)
      call Var_Mod_Allocate_Gradients(f22)

      if(turbulence_statistics) then
        call Var_Mod_Allocate_Statistics(f22)
      end if
    end if

    ! Other variables such as time scale, length scale and production
    allocate(t_scale(-grid % n_bnd_cells:grid % n_cells));  t_scale = 0.
    allocate(l_scale(-grid % n_bnd_cells:grid % n_cells));  l_scale = 0.
    allocate(p_kin  (-grid % n_bnd_cells:grid % n_cells));  p_kin   = 0.
    if(turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
      allocate(eps_tot(-grid % n_bnd_cells:grid % n_cells)); eps_tot = 0.
    end if

  end if

  !----------------------!
  !   Spalart Allmaras   !
  !----------------------!
  if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
     turbulence_model .eq. DES_SPALART) then
    call Var_Mod_Allocate_Solution('VIS', vis, grid)

    if(turbulence_statistics) then
      call Var_Mod_Allocate_Statistics(vis)
    end if
  end if

  !-----------------------!
  !   Smagorinsky model   !
  !-----------------------!
  if(turbulence_model .eq. LES_SMAGORINSKY) then
    allocate(nearest_wall_cell(-grid % n_bnd_cells:grid % n_cells))
    nearest_wall_cell = 0
  end if

  !----------------!
  !   Wale model   !
  !----------------!
  if(turbulence_model .eq. LES_WALE) then
    allocate(wale_v(-grid % n_bnd_cells:grid % n_cells));  wale_v = 0.
  end if

  !-------------------!
  !   Dynamic model   !
  !-------------------!
  if(turbulence_model .eq. LES_DYNAMIC) then
    allocate(c_dyn(-grid % n_bnd_cells:grid % n_cells));  c_dyn = 0.
  end if

  !-----------------------------------------------------------------------!
  !                                                                       !
  !   Turbulent statistics for all models except second moment closures   !
  !                                                                       !
  !-----------------------------------------------------------------------!
  if(turbulence_statistics) then

    ! For second moment closures, memory for statistics was allocated above
    if(turbulence_model .ne. RSM_HANJALIC_JAKIRLIC .and.  &
       turbulence_model .ne. RSM_MANCEAU_HANJALIC) then

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

    end if

  end if

  !-----------------------------!
  !                             !
  !   User scalars and arrays   !
  !                             !
  !-----------------------------!
  call Control_Mod_Number_Of_User_Scalars(n_user_scalars, verbose = .true.)
  call Control_Mod_Number_Of_User_Arrays (n_user_arrays,  verbose = .true.)
  call User_Mod_Allocate(grid)

  end subroutine

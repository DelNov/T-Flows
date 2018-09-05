!==============================================================================!
  subroutine Backup_Mod_Save(grid, time_step, name_save)
!------------------------------------------------------------------------------!
!   Saves backup files name.backup                                             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Name_Mod, only: problem_name
  use Const_Mod
  use Flow_Mod
  use Comm_Mod
  use Rans_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  integer          :: time_step
  character(len=*) :: name_save
!-----------------------------------[Locals]-----------------------------------!
  character(len=80)    :: name_out, store_name
  integer              :: fh, d
  integer              :: s, c, c1, c2, cm1, cm2, mc, max_cnt
  integer, allocatable :: cell_map(:)
  integer, allocatable :: face_cnt(:)
  integer, allocatable :: face_pnt(:,:)
  real,    allocatable :: cell_flx(:)
  character(len=12)    :: cf_name = 'cell_flux_00'
!==============================================================================!

  store_name = problem_name

  problem_name = name_save

  ! Name backup file
  call Name_File(0, name_out, '.backup')

  ! Open backup file
  call Comm_Mod_Open_File_Write(fh, name_out)

  ! Create new types
  call Comm_Mod_Create_New_Types(grid)

  ! Initialize displacement
  d = 0

  !-----------------------------------------------------------------------!
  !   Save cell-centre coordinates.  Could be useful for interpolations   !
  !-----------------------------------------------------------------------!
  call Backup_Mod_Write_Cell_Bnd(fh, d, 'x_coordinates', grid % xc(-nb_s:nc_s))
  call Backup_Mod_Write_Cell_Bnd(fh, d, 'y_coordinates', grid % yc(-nb_s:nc_s))
  call Backup_Mod_Write_Cell_Bnd(fh, d, 'z_coordinates', grid % zc(-nb_s:nc_s))

  !---------------!
  !               !
  !   Save data   !
  !               !
  !---------------!

  ! Time step
  call Backup_Mod_Write_Int(fh, d, 'time_step', time_step)

  ! Number of processors
  call Backup_Mod_Write_Int(fh, d, 'n_proc', n_proc)

  ! Bulk flows and pressure drops in each direction
  call Backup_Mod_Write_Real(fh, d, 'bulk_flux_x',   bulk % flux_x)
  call Backup_Mod_Write_Real(fh, d, 'bulk_flux_y',   bulk % flux_y)
  call Backup_Mod_Write_Real(fh, d, 'bulk_flux_z',   bulk % flux_z)
  call Backup_Mod_Write_Real(fh, d, 'bulk_p_drop_x', bulk % p_drop_x)
  call Backup_Mod_Write_Real(fh, d, 'bulk_p_drop_y', bulk % p_drop_y)
  call Backup_Mod_Write_Real(fh, d, 'bulk_p_drop_z', bulk % p_drop_z)

  !----------------------------!
  !                            !
  !   Navier-Stokes equation   !
  !                            !
  !----------------------------!

  !--------------!
  !   Velocity   !
  !--------------!
  call Backup_Mod_Write_Variable(fh, d, 'u_velocity', u)
  call Backup_Mod_Write_Variable(fh, d, 'v_velocity', v)
  call Backup_Mod_Write_Variable(fh, d, 'w_velocity', w)

  !--------------------------------------!
  !   Pressure and pressure correction   !
  !--------------------------------------!
  call Backup_Mod_Write_Cell_Bnd(fh, d, 'press',      p % n(-nb_s:nc_s))
  call Backup_Mod_Write_Cell_Bnd(fh, d, 'press_corr',pp % n(-nb_s:nc_s))

  !----------------------!
  !   Mass flow raters   !
  !----------------------!
  call Backup_Mod_Write_Face(grid, flux, fh, d) 

  !--------------!
  !              !
  !   Etnhalpy   !
  !              !
  !--------------!
  if(heat_transfer) then
    call Backup_Mod_Write_Variable(fh, d, 'temp', t)
  end if

  !-----------------------!
  !                       !
  !   Turbulence models   !
  !                       !
  !-----------------------!

  !-----------------!
  !   K-eps model   !
  !-----------------!
  if(turbulence_model .eq. K_EPS) then

    ! K and epsilon
    call Backup_Mod_Write_Variable(fh, d, 'kin', kin)
    call Backup_Mod_Write_Variable(fh, d, 'eps', eps)

    ! Other turbulent quantities
    call Backup_Mod_Write_Cell_Bnd(fh, d, 'p_kin',    p_kin   (-nb_s:nc_s))
    call Backup_Mod_Write_Cell_Bnd(fh, d, 'u_tau',    u_tau   (-nb_s:nc_s))
    call Backup_Mod_Write_Cell_Bnd(fh, d, 'y_plus',   y_plus  (-nb_s:nc_s))
    call Backup_Mod_Write_Cell_Bnd(fh, d, 'vis_t',    vis_t   (-nb_s:nc_s))
    call Backup_Mod_Write_Cell_Bnd(fh, d, 'vis_wall', vis_wall(-nb_s:nc_s))
    call Backup_Mod_Write_Cell    (fh, d, 'tau_wall', tau_wall  (1:nc_s))
  end if

  !------------------------!
  !   K-eps-zeta-f model   !
  !------------------------!
  if(turbulence_model .eq. K_EPS_ZETA_F) then

    ! K, eps, zeta and f22
    call Backup_Mod_Write_Variable(fh, d, 'kin',  kin)
    call Backup_Mod_Write_Variable(fh, d, 'eps',  eps)
    call Backup_Mod_Write_Variable(fh, d, 'zeta', zeta)
    call Backup_Mod_Write_Variable(fh, d, 'f22',  f22)

    ! Other turbulent quantities
    call Backup_Mod_Write_Cell_Bnd(fh, d, 'p_kin',    p_kin   (-nb_s:nc_s))
    call Backup_Mod_Write_Cell_Bnd(fh, d, 'u_tau',    u_tau   (-nb_s:nc_s))
    call Backup_Mod_Write_Cell_Bnd(fh, d, 'y_plus',   y_plus  (-nb_s:nc_s))
    call Backup_Mod_Write_Cell_Bnd(fh, d, 'vis_t',    vis_t   (-nb_s:nc_s))
    call Backup_Mod_Write_Cell_Bnd(fh, d, 'vis_wall', vis_wall(-nb_s:nc_s))
    call Backup_Mod_Write_Cell    (fh, d, 'tau_wall', tau_wall  (1:nc_s))
    call Backup_Mod_Write_Cell_Bnd(fh, d, 't_scale',  t_scale(-nb_s:nc_s))
    call Backup_Mod_Write_Cell_Bnd(fh, d, 'l_scale',  l_scale(-nb_s:nc_s))
  end if

  !----------------------------!
  !   Reynolds stress models   !
  !----------------------------!
  if(turbulence_model .eq. REYNOLDS_STRESS .or.  &
     turbulence_model .eq. HANJALIC_JAKIRLIC) then

    ! Reynolds stresses
    call Backup_Mod_Write_Variable(fh, d, 'uu',  uu)
    call Backup_Mod_Write_Variable(fh, d, 'vv',  vv)
    call Backup_Mod_Write_Variable(fh, d, 'ww',  ww)
    call Backup_Mod_Write_Variable(fh, d, 'uv',  uv)
    call Backup_Mod_Write_Variable(fh, d, 'uw',  uw)
    call Backup_Mod_Write_Variable(fh, d, 'vw',  vw)

    ! Epsilon
    call Backup_Mod_Write_Variable(fh, d, 'eps', eps)

    ! F22
    if(turbulence_model .eq. HANJALIC_JAKIRLIC) then
      call Backup_Mod_Write_Variable(fh, d, 'f22',  f22)
    end if

    ! Other turbulent quantities 
    call Backup_Mod_Write_Cell_Bnd(fh, d, 'vis_t', vis_t(-nb_s:nc_s))
  end if

  !-----------------------------------------!
  !                                         !
  !   Turbulent statistics for all models   !
  !                                         !
  !-----------------------------------------!
  if(turbulence_statistics) then
    call Backup_Mod_Write_Variable_Mean(fh, d, 'u_mean', u)
    call Backup_Mod_Write_Variable_Mean(fh, d, 'v_mean', v)
    call Backup_Mod_Write_Variable_Mean(fh, d, 'w_mean', w)

    call Backup_Mod_Write_Variable_Mean(fh, d, 'uu_mean', uu)
    call Backup_Mod_Write_Variable_Mean(fh, d, 'vv_mean', vv)
    call Backup_Mod_Write_Variable_Mean(fh, d, 'ww_mean', ww)
    call Backup_Mod_Write_Variable_Mean(fh, d, 'uv_mean', uv)
    call Backup_Mod_Write_Variable_Mean(fh, d, 'uw_mean', uw)
    call Backup_Mod_Write_Variable_Mean(fh, d, 'vw_mean', vw)

    if(heat_transfer) then
      call Backup_Mod_Write_Variable_Mean(fh, d, 't_mean',  t)
      call Backup_Mod_Write_Variable_Mean(fh, d, 'tt_mean', tt)
      call Backup_Mod_Write_Variable_Mean(fh, d, 'ut_mean', ut)
      call Backup_Mod_Write_Variable_Mean(fh, d, 'vt_mean', vt)
      call Backup_Mod_Write_Variable_Mean(fh, d, 'wt_mean', wt)
    end if
  end if

  !------------------------------!
  !                              !
  !   User scalars are missing   !
  !                              !
  !------------------------------!

  ! Close backup file
  call Comm_Mod_Close_File(fh)

  problem_name = store_name

  end subroutine

!==============================================================================!
  subroutine Backup_Mod_Load(fld, swr, time_step, time_step_stat, backup)
!------------------------------------------------------------------------------!
!   Loads backup files name.backup                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: fld
  type(Swarm_Type), target :: swr
  integer                  :: time_step       ! current time step
  integer                  :: time_step_stat  ! starting step for statistics
  logical                  :: backup, present
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Bulk_Type), pointer :: bulk
  character(len=80)        :: name_in, answer
  integer                  :: fh, d, vc
!==============================================================================!

  ! Take aliases
  grid => fld % pnt_grid
  bulk => fld % bulk

  ! Full name is specified in control file
  call Control_Mod_Load_Backup_Name(name_in)

  answer = name_in
  call To_Upper_Case(answer)

  backup = .true.
  if(answer .eq. 'SKIP') then
    backup = .false.
    return
  end if

  inquire(file=trim(name_in), exist=present )
  if(.not.present) then
    if(this_proc < 2) then
      print *, "# ERROR!  Backup file ", trim(name_in), " was not found."
      print *, "# Exiting!"
    end if
    call Comm_Mod_End
  end if

  ! Open backup file
  call Comm_Mod_Open_File_Read(fh, name_in)

  ! Create new types
  call Comm_Mod_Create_New_Types(grid)

  ! Initialize displacement and variable count
  d  = 0
  vc = 0

  !---------------------------------------------!
  !   Variable count - important for checking   !
  !---------------------------------------------!
  call Backup_Mod_Read_Int(fh, d, 2048, 'variable_count', vc)
  if(vc .eq. 0) vc = 2048  ! for backward compatibility

  if(this_proc < 2) then
    print *, "# Backup file holds ", vc, " variables."
  end if

  !---------------!
  !               !
  !   Load data   !
  !               !
  !---------------!

  !-----------------------------------------------!
  !   Skip three coordinates for the time being   !
  !-----------------------------------------------!
  ! call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'x_coords', grid % xc(-nb_s:nc_s))
  ! call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'y_coords', grid % yc(-nb_s:nc_s))
  ! call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'z_coords', grid % zc(-nb_s:nc_s))

  ! Time step
  call Backup_Mod_Read_Int(fh, d, vc, 'time_step', time_step)

  ! Bulk flows and pressure drops in each direction
  call Backup_Mod_Read_Real(fh, d, vc, 'bulk_flux_x',   bulk % flux_x)
  call Backup_Mod_Read_Real(fh, d, vc, 'bulk_flux_y',   bulk % flux_y)
  call Backup_Mod_Read_Real(fh, d, vc, 'bulk_flux_z',   bulk % flux_z)
  call Backup_Mod_Read_Real(fh, d, vc, 'bulk_flux_x_o', bulk % flux_x_o)
  call Backup_Mod_Read_Real(fh, d, vc, 'bulk_flux_y_o', bulk % flux_y_o)
  call Backup_Mod_Read_Real(fh, d, vc, 'bulk_flux_z_o', bulk % flux_z_o)
  call Backup_Mod_Read_Real(fh, d, vc, 'bulk_p_drop_x', bulk % p_drop_x)
  call Backup_Mod_Read_Real(fh, d, vc, 'bulk_p_drop_y', bulk % p_drop_y)
  call Backup_Mod_Read_Real(fh, d, vc, 'bulk_p_drop_z', bulk % p_drop_z)

  !----------------------------!
  !                            !
  !   Navier-Stokes equation   !
  !                            !
  !----------------------------!

  !--------------!
  !   Velocity   !
  !--------------!
  call Backup_Mod_Read_Variable(fh, d, vc, 'u_velocity', fld % u)
  call Backup_Mod_Read_Variable(fh, d, vc, 'v_velocity', fld % v)
  call Backup_Mod_Read_Variable(fh, d, vc, 'w_velocity', fld % w)

  !--------------------------------------!
  !   Pressure and pressure correction   !
  !--------------------------------------!
  call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'press',      fld % p % n(-nb_s:nc_s))
  call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'press_corr', fld %pp % n(-nb_s:nc_s))

  !----------------------!
  !   Mass flow raters   !
  !----------------------!
  call Backup_Mod_Read_Face(fh, d, vc, grid, fld % flux)

  !--------------!
  !              !
  !   Etnhalpy   !
  !              !
  !--------------!
  if(heat_transfer) then
    call Backup_Mod_Read_Variable(fh, d, vc, 'temp', fld % t)
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
    call Backup_Mod_Read_Variable(fh, d, vc, 'kin', kin)
    call Backup_Mod_Read_Variable(fh, d, vc, 'eps', eps)

    ! Other turbulent quantities
    call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'p_kin',    p_kin   (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'u_tau',    u_tau   (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'y_plus',   y_plus  (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'vis_t',    vis_t   (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'vis_wall', vis_wall(-nb_s:nc_s))
    call Backup_Mod_Read_Cell    (fh, d, vc, 'tau_wall', tau_wall  (1:nc_s)  )

    ! Turbulence quantities connected with heat transfer
    if(heat_transfer) then
      call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'con_wall', con_wall(-nb_s:nc_s))
    end if
  end if

  !------------------------!
  !   K-eps-zeta-f model   !
  !------------------------!
  if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
     turbulence_model .eq. HYBRID_LES_RANS) then

    ! K, eps, zeta and f22
    call Backup_Mod_Read_Variable(fh, d, vc, 'kin',  kin)
    call Backup_Mod_Read_Variable(fh, d, vc, 'eps',  eps)
    call Backup_Mod_Read_Variable(fh, d, vc, 'zeta', zeta)
    call Backup_Mod_Read_Variable(fh, d, vc, 'f22',  f22)

    ! Other turbulent quantities
    call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'p_kin',    p_kin   (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'u_tau',    u_tau   (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'y_plus',   y_plus  (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'vis_t',    vis_t   (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'vis_wall', vis_wall(-nb_s:nc_s))
    call Backup_Mod_Read_Cell    (fh, d, vc, 'tau_wall', tau_wall  (1:nc_s)  )
    call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 't_scale',  t_scale (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'l_scale',  l_scale (-nb_s:nc_s))

    ! Turbulence quantities connected with heat transfer
  end if

  if(turbulence_model .eq. K_EPS_ZETA_F .and. heat_transfer) then
    call Backup_Mod_Read_Variable(fh, d, vc, 't2',       t2)
    call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'p_t2',     p_t2    (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'con_wall', con_wall(-nb_s:nc_s))
  else if (turbulence_model .eq. HYBRID_LES_RANS .and. heat_transfer) then
    call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'con_wall', con_wall(-nb_s:nc_s))
  end if

  !----------------------------!
  !   Reynolds stress models   !
  !----------------------------!
  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then

    ! Reynolds stresses
    call Backup_Mod_Read_Variable(fh, d, vc, 'uu',  uu)
    call Backup_Mod_Read_Variable(fh, d, vc, 'vv',  vv)
    call Backup_Mod_Read_Variable(fh, d, vc, 'ww',  ww)
    call Backup_Mod_Read_Variable(fh, d, vc, 'uv',  uv)
    call Backup_Mod_Read_Variable(fh, d, vc, 'uw',  uw)
    call Backup_Mod_Read_Variable(fh, d, vc, 'vw',  vw)

    ! Epsilon
    call Backup_Mod_Read_Variable(fh, d, vc, 'eps', eps)

    ! F22
    if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
      call Backup_Mod_Read_Variable(fh, d, vc, 'f22',  f22)
    end if

    ! Other turbulent quantities ?
    call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'vis_t', vis_t(-nb_s:nc_s))

    ! Turbulence quantities connected with heat transfer
    if(heat_transfer) then
      call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'con_wall', con_wall(-nb_s:nc_s))
    end if
  end if

  !-----------------------------------------!
  !                                         !
  !   Turbulent statistics for all models   !
  !                                         !
  !-----------------------------------------!
  if(turbulence_statistics .and.  &
     time_step > time_step_stat) then

    call Backup_Mod_Read_Variable_Mean(fh, d, vc, 'u_mean', fld % u)
    call Backup_Mod_Read_Variable_Mean(fh, d, vc, 'v_mean', fld % v)
    call Backup_Mod_Read_Variable_Mean(fh, d, vc, 'w_mean', fld % w)

    call Backup_Mod_Read_Variable_Mean(fh, d, vc, 'uu_mean', uu)
    call Backup_Mod_Read_Variable_Mean(fh, d, vc, 'vv_mean', vv)
    call Backup_Mod_Read_Variable_Mean(fh, d, vc, 'ww_mean', ww)
    call Backup_Mod_Read_Variable_Mean(fh, d, vc, 'uv_mean', uv)
    call Backup_Mod_Read_Variable_Mean(fh, d, vc, 'uw_mean', uw)
    call Backup_Mod_Read_Variable_Mean(fh, d, vc, 'vw_mean', vw)

    if(heat_transfer) then
      call Backup_Mod_Read_Variable_Mean(fh, d, vc, 't_mean',  fld % t)
      call Backup_Mod_Read_Variable_Mean(fh, d, vc, 'tt_mean', tt)
      call Backup_Mod_Read_Variable_Mean(fh, d, vc, 'ut_mean', ut)
      call Backup_Mod_Read_Variable_Mean(fh, d, vc, 'vt_mean', vt)
      call Backup_Mod_Read_Variable_Mean(fh, d, vc, 'wt_mean', wt)
    end if
  end if

  !--------------------------!
  !                          !
  !   Swarm (of particles)   !
  !                          !
  !--------------------------!
  call Backup_Mod_Read_Swarm(fh, d, vc, swr)

  !------------------------------!
  !                              !
  !   User scalars are missing   !
  !                              !
  !------------------------------!

  ! Close backup file
  call Comm_Mod_Close_File(fh)

  end subroutine

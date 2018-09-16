!==============================================================================!
  subroutine Backup_Mod_Load(grid, time_step, backup)
!------------------------------------------------------------------------------!
!   Loads backup files name.backup                                             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
! use Les_Mod
  use Comm_Mod
  use Rans_Mod
  use Grid_Mod
  use Grad_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: time_step
  logical         :: backup, present
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: name_in, answer
  integer           :: fh, d
!==============================================================================!

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
      print *, "# Backup file ", trim(name_in)," was not found.  Exiting!"
    end if
    stop
  end if

  ! Open backup file
  call Comm_Mod_Open_File_Read(fh, name_in)

  ! Create new types
  call Comm_Mod_Create_New_Types(grid)

  ! Initialize displacement
  d = 0

  !-----------------------------------------------!
  !   Skip three coordinates for the time being   !
  !-----------------------------------------------!
  ! call Backup_Mod_Read_Cell_Bnd(fh, d, 'x_coordinates', grid % xc(-nb_s:nc_s))
  ! call Backup_Mod_Read_Cell_Bnd(fh, d, 'y_coordinates', grid % yc(-nb_s:nc_s))
  ! call Backup_Mod_Read_Cell_Bnd(fh, d, 'z_coordinates', grid % zc(-nb_s:nc_s))

  !---------------!
  !               !
  !   Load data   !
  !               !
  !---------------!

  ! Time step
  call Backup_Mod_Read_Int(fh, d, 'time_step', time_step)

  ! Bulk flows and pressure drops in each direction
  call Backup_Mod_Read_Real(fh, d, 'bulk_flux_x',   bulk % flux_x)
  call Backup_Mod_Read_Real(fh, d, 'bulk_flux_y',   bulk % flux_y)
  call Backup_Mod_Read_Real(fh, d, 'bulk_flux_z',   bulk % flux_z)
  call Backup_Mod_Read_Real(fh, d, 'bulk_p_drop_x', bulk % p_drop_x)
  call Backup_Mod_Read_Real(fh, d, 'bulk_p_drop_y', bulk % p_drop_y)
  call Backup_Mod_Read_Real(fh, d, 'bulk_p_drop_z', bulk % p_drop_z)

  !----------------------------!
  !                            !
  !   Navier-Stokes equation   !
  !                            !
  !----------------------------!

  !--------------!
  !   Velocity   !
  !--------------!
  call Backup_Mod_Read_Variable(fh, d, 'u_velocity', u)
  call Backup_Mod_Read_Variable(fh, d, 'v_velocity', v)
  call Backup_Mod_Read_Variable(fh, d, 'w_velocity', w)

  !--------------------------------------!
  !   Pressure and pressure correction   !
  !--------------------------------------!
  call Backup_Mod_Read_Cell_Bnd(fh, d, 'press',       p % n(-nb_s:nc_s))
  call Backup_Mod_Read_Cell_Bnd(fh, d, 'press_corr', pp % n(-nb_s:nc_s))

  !----------------------!
  !   Mass flow raters   !
  !----------------------!
  call Backup_Mod_Read_Face(grid, flux, fh, d) 

  !--------------!
  !              !
  !   Etnhalpy   !
  !              !
  !--------------!
  if(heat_transfer) then
    call Backup_Mod_Read_Variable(fh, d, 'temp', t)
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
    call Backup_Mod_Read_Variable(fh, d, 'kin', kin)
    call Backup_Mod_Read_Variable(fh, d, 'eps', eps)

    ! Other turbulent quantities
    call Backup_Mod_Read_Cell_Bnd(fh, d, 'p_kin',    p_kin   (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh, d, 'u_tau',    u_tau   (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh, d, 'y_plus',   y_plus  (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh, d, 'vis_t',    vis_t   (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh, d, 'vis_wall', vis_wall(-nb_s:nc_s))
    call Backup_Mod_Read_Cell    (fh, d, 'tau_wall', tau_wall  (1:nc_s)  )
  end if

  !------------------------!
  !   K-eps-zeta-f model   !
  !------------------------!
  if(turbulence_model .eq. K_EPS_ZETA_F) then

    ! K, eps, zeta and f22
    call Backup_Mod_Read_Variable(fh, d, 'kin',  kin)
    call Backup_Mod_Read_Variable(fh, d, 'eps',  eps)
    call Backup_Mod_Read_Variable(fh, d, 'zeta', zeta)
    call Backup_Mod_Read_Variable(fh, d, 'f22',  f22)

    ! Other turbulent quantities
    call Backup_Mod_Read_Cell_Bnd(fh, d, 'p_kin',    p_kin   (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh, d, 'u_tau',    u_tau   (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh, d, 'y_plus',   y_plus  (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh, d, 'vis_t',    vis_t   (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh, d, 'vis_wall', vis_wall(-nb_s:nc_s))
    call Backup_Mod_Read_Cell    (fh, d, 'tau_wall', tau_wall  (1:nc_s)  )
    call Backup_Mod_Read_Cell_Bnd(fh, d, 't_scale',  t_scale (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh, d, 'l_scale',  l_scale (-nb_s:nc_s))
  end if

  !----------------------------!
  !   Reynolds stress models   !
  !----------------------------!
  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then

    ! Reynolds stresses
    call Backup_Mod_Read_Variable(fh, d, 'uu',  uu)
    call Backup_Mod_Read_Variable(fh, d, 'vv',  vv)
    call Backup_Mod_Read_Variable(fh, d, 'ww',  ww)
    call Backup_Mod_Read_Variable(fh, d, 'uv',  uv)
    call Backup_Mod_Read_Variable(fh, d, 'uw',  uw)
    call Backup_Mod_Read_Variable(fh, d, 'vw',  vw)

    ! Epsilon
    call Backup_Mod_Read_Variable(fh, d, 'eps', eps)

    ! F22
    if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
      call Backup_Mod_Read_Variable(fh, d, 'f22',  f22)
    end if

    ! Other turbulent quantities ?
    call Backup_Mod_Read_Cell_Bnd(fh, d, 'vis_t', vis_t(-nb_s:nc_s))
  end if

  !-----------------------------------------!
  !                                         !
  !   Turbulent statistics for all models   !
  !                                         !
  !-----------------------------------------!
  if(turbulence_statistics) then
    call Backup_Mod_Read_Variable_Mean(fh, d, 'u_mean', u)
    call Backup_Mod_Read_Variable_Mean(fh, d, 'v_mean', v)
    call Backup_Mod_Read_Variable_Mean(fh, d, 'w_mean', w)

    call Backup_Mod_Read_Variable_Mean(fh, d, 'uu_mean', uu)
    call Backup_Mod_Read_Variable_Mean(fh, d, 'vv_mean', vv)
    call Backup_Mod_Read_Variable_Mean(fh, d, 'ww_mean', ww)
    call Backup_Mod_Read_Variable_Mean(fh, d, 'uv_mean', uv)
    call Backup_Mod_Read_Variable_Mean(fh, d, 'uw_mean', uw)
    call Backup_Mod_Read_Variable_Mean(fh, d, 'vw_mean', vw)

    if(heat_transfer) then
      call Backup_Mod_Read_Variable_Mean(fh, d, 't_mean',  t)
      call Backup_Mod_Read_Variable_Mean(fh, d, 'tt_mean', tt)
      call Backup_Mod_Read_Variable_Mean(fh, d, 'ut_mean', ut)
      call Backup_Mod_Read_Variable_Mean(fh, d, 'vt_mean', vt)
      call Backup_Mod_Read_Variable_Mean(fh, d, 'wt_mean', wt)
    end if
  end if

  !------------------------------!
  !                              !
  !   User scalars are missing   !
  !                              !
  !------------------------------!

  ! Close backup file
  call Comm_Mod_Close_File(fh)

  end subroutine

!==============================================================================!
  subroutine Load_Backup(grid, time_step, restart)
!------------------------------------------------------------------------------!
!   Read backup files. name.backup                                             !
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
  logical         :: restart, present
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: name_in, answer
  integer           :: fh, d
!==============================================================================!

  ! Full name is specified in control file
  call Control_Mod_Load_Backup_Name(name_in)

  answer = name_in
  call To_Upper_Case(answer)
  if(answer .eq. 'SKIP') then
    restart = .false.
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
  call Comm_Mod_Create_New_Types()

  ! Initialize displacement
  d = 0

  !-----------------------------------------------!
  !   Skip three coordinates for the time being   !
  !-----------------------------------------------!
  ! call Read_Backup_3_Cell_Bnd(fh, d, 'coordinates',   &
  !                             grid % xc(-nb_s:nc_s),  &
  !                             grid % yc(-nb_s:nc_s),  &
  !                             grid % zc(-nb_s:nc_s))

  !---------------!
  !               !
  !   Load data   !
  !               !
  !---------------!

  ! Time step
  call Read_Backup_Int(fh, d, 'time_step', time_step)

  ! Bulk flows and pressure drops in each direction
  call Read_Backup_Real(fh, d, 'bulk_flux_x',   bulk(1) % flux_x)
  call Read_Backup_Real(fh, d, 'bulk_flux_y',   bulk(1) % flux_y)
  call Read_Backup_Real(fh, d, 'bulk_flux_z',   bulk(1) % flux_z)
  call Read_Backup_Real(fh, d, 'bulk_p_drop_x', bulk(1) % p_drop_x)
  call Read_Backup_Real(fh, d, 'bulk_p_drop_y', bulk(1) % p_drop_y)
  call Read_Backup_Real(fh, d, 'bulk_p_drop_z', bulk(1) % p_drop_z)

  !----------------------------!
  !                            !
  !   Navier-Stokes equation   !
  !                            !
  !----------------------------!

  !--------------!
  !   Velocity   !
  !--------------!
  call Read_Backup_Variable(fh, d, 'u_velocity', u)
  call Read_Backup_Variable(fh, d, 'v_velocity', v)
  call Read_Backup_Variable(fh, d, 'w_velocity', w)

  !--------------------------------------!
  !   Pressure and pressure correction   !
  !--------------------------------------!
  call Read_Backup_Cell_Bnd(fh, d, 'press',       p  % n(-nb_s:nc_s))
  call Read_Backup_Cell_Bnd(fh, d, 'press_corr',  pp % n(-nb_s:nc_s))

  !----------------------------------------------------!
  !   Mass flow rates (ask Egor if name is correct?)   !
  !----------------------------------------------------!
  call Read_Backup_Face(fh, d, 'mass_flow_rate', flux(1:nf_s+nbf_s))
  call Calculate_Mass_Flow_Rate(grid)

  !--------------!
  !              !
  !   Etnhalpy   !
  !              !
  !--------------!
  if(heat_transfer .eq. YES) then
    call Read_Backup_Variable(fh, d, 'temp', t)
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
    call Read_Backup_Variable(fh, d, 'kin', kin)
    call Read_Backup_Variable(fh, d, 'eps', eps)

    ! Other turbulent quantities
    call Read_Backup_Cell_Bnd(fh, d, 'p_kin',    p_kin   (-nb_s:nc_s))
    call Read_Backup_Cell_Bnd(fh, d, 'u_tau',    u_tau   (-nb_s:nc_s))
    call Read_Backup_Cell_Bnd(fh, d, 'y_plus',   y_plus  (-nb_s:nc_s))
    call Read_Backup_Cell_Bnd(fh, d, 'vis_wall', vis_wall(-nb_s:nc_s))
    call Read_Backup_Cell    (fh, d, 'tau_wall', tau_wall  (1:nc_s)  )
  end if

  !------------------------!
  !   K-eps-zeta-f model   !
  !------------------------!
  if(turbulence_model .eq. K_EPS_ZETA_F) then

    ! K, eps, zeta and f22
    call Read_Backup_Variable(fh, d, 'kin',  kin)
    call Read_Backup_Variable(fh, d, 'eps',  eps)
    call Read_Backup_Variable(fh, d, 'zeta', zeta)
    call Read_Backup_Variable(fh, d, 'f22',  f22)

    ! Other turbulent quantities
    call Read_Backup_Cell_Bnd(fh, d, 'p_kin',    p_kin   (-nb_s:nc_s))
    call Read_Backup_Cell_Bnd(fh, d, 'u_tau',    u_tau   (-nb_s:nc_s))
    call Read_Backup_Cell_Bnd(fh, d, 'y_plus',   y_plus  (-nb_s:nc_s))
    call Read_Backup_Cell_Bnd(fh, d, 'vis_wall', vis_wall(-nb_s:nc_s))
    call Read_Backup_Cell    (fh, d, 'tau_wall', tau_wall  (1:nc_s)  )
    call Read_Backup_Cell_Bnd(fh, d, 't_scale',  t_scale (-nb_s:nc_s))
    call Read_Backup_Cell_Bnd(fh, d, 'l_scale',  l_scale (-nb_s:nc_s))
  end if

  !----------------------------!
  !   Reynolds stress models   !
  !----------------------------!
  if(turbulence_model .eq. REYNOLDS_STRESS .or.  &
     turbulence_model .eq. HANJALIC_JAKIRLIC) then

    ! Reynolds stresses
    call Read_Backup_Variable(fh, d, 'uu',  uu)
    call Read_Backup_Variable(fh, d, 'vv',  vv)
    call Read_Backup_Variable(fh, d, 'ww',  ww)
    call Read_Backup_Variable(fh, d, 'uv',  uv)
    call Read_Backup_Variable(fh, d, 'uw',  uw)
    call Read_Backup_Variable(fh, d, 'vw',  vw)

    ! Epsilon
    call Read_Backup_Variable(fh, d, 'eps', eps)

    ! F22
    if(turbulence_model .eq. REYNOLDS_STRESS) then
      call Read_Backup_Variable(fh, d, 'f22',  f22)
    end if

    ! Other turbulent quantities ?
  end if

  !-----------------------------------------!
  !                                         !
  !   Turbulent statistics for all models   !
  !                                         !
  !-----------------------------------------!
  if(turbulence_statistics .eq. YES) then
    call Read_Backup_Variable_Mean(fh, d, 'u_mean', u)
    call Read_Backup_Variable_Mean(fh, d, 'v_mean', v)
    call Read_Backup_Variable_Mean(fh, d, 'w_mean', w)

    call Read_Backup_Variable_Mean(fh, d, 'uu_mean', uu)
    call Read_Backup_Variable_Mean(fh, d, 'vv_mean', vv)
    call Read_Backup_Variable_Mean(fh, d, 'ww_mean', ww)
    call Read_Backup_Variable_Mean(fh, d, 'uv_mean', uv)
    call Read_Backup_Variable_Mean(fh, d, 'uw_mean', uw)
    call Read_Backup_Variable_Mean(fh, d, 'vw_mean', vw)

    if(heat_transfer .eq. YES) then
      call Read_Backup_Variable_Mean(fh, d, 't_mean',  t)
      call Read_Backup_Variable_Mean(fh, d, 'tt_mean', tt)
      call Read_Backup_Variable_Mean(fh, d, 'ut_mean', ut)
      call Read_Backup_Variable_Mean(fh, d, 'vt_mean', vt)
      call Read_Backup_Variable_Mean(fh, d, 'wt_mean', wt)
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
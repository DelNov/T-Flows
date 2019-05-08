!==============================================================================!
  subroutine Backup_Mod_Load(fld, swr, tur,  &
                             time_step, time_step_stat, backup)
!------------------------------------------------------------------------------!
!   Loads backup files name.backup                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: fld
  type(Swarm_Type), target :: swr
  type(Turb_Type),  target :: tur
  integer                  :: time_step       ! current time step
  integer                  :: time_step_stat  ! starting step for statistics
  logical                  :: backup, present
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Bulk_Type), pointer :: bulk
  character(len=80)        :: name_in, answer
  integer                  :: fh,d,vc
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
  ! call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'x_coords', grid % xc(-nb_s:nc_s))
  ! call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'y_coords', grid % yc(-nb_s:nc_s))
  ! call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'z_coords', grid % zc(-nb_s:nc_s))

  ! Time step
  call Backup_Mod_Read_Int(fh,d,vc, 'time_step', time_step)

  ! Bulk flows and pressure drops in each direction
  call Backup_Mod_Read_Real(fh,d,vc, 'bulk_flux_x',   bulk % flux_x)
  call Backup_Mod_Read_Real(fh,d,vc, 'bulk_flux_y',   bulk % flux_y)
  call Backup_Mod_Read_Real(fh,d,vc, 'bulk_flux_z',   bulk % flux_z)
  call Backup_Mod_Read_Real(fh,d,vc, 'bulk_flux_x_o', bulk % flux_x_o)
  call Backup_Mod_Read_Real(fh,d,vc, 'bulk_flux_y_o', bulk % flux_y_o)
  call Backup_Mod_Read_Real(fh,d,vc, 'bulk_flux_z_o', bulk % flux_z_o)
  call Backup_Mod_Read_Real(fh,d,vc, 'bulk_p_drop_x', bulk % p_drop_x)
  call Backup_Mod_Read_Real(fh,d,vc, 'bulk_p_drop_y', bulk % p_drop_y)
  call Backup_Mod_Read_Real(fh,d,vc, 'bulk_p_drop_z', bulk % p_drop_z)

  !----------------------------!
  !                            !
  !   Navier-Stokes equation   !
  !                            !
  !----------------------------!

  !--------------!
  !   Velocity   !
  !--------------!
  call Backup_Mod_Read_Variable(fh,d,vc, 'u_velocity', fld % u)
  call Backup_Mod_Read_Variable(fh,d,vc, 'v_velocity', fld % v)
  call Backup_Mod_Read_Variable(fh,d,vc, 'w_velocity', fld % w)

  !--------------------------------------!
  !   Pressure and pressure correction   !
  !--------------------------------------!
  call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'press',      fld % p % n(-nb_s:nc_s))
  call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'press_corr', fld %pp % n(-nb_s:nc_s))

  !----------------------!
  !   Mass flow raters   !
  !----------------------!
  call Backup_Mod_Read_Face(fh,d,vc, grid, fld % flux)

  !--------------!
  !              !
  !   Etnhalpy   !
  !              !
  !--------------!
  if(heat_transfer) then
    call Backup_Mod_Read_Variable(fh,d,vc, 'temp', fld % t)
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
    call Backup_Mod_Read_Variable(fh,d,vc, 'kin', tur % kin)
    call Backup_Mod_Read_Variable(fh,d,vc, 'eps', tur % eps)

    ! Other turbulent quantities
    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'p_kin',    p_kin   (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'y_plus',   y_plus  (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'vis_t',    vis_t   (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'vis_wall', vis_wall(-nb_s:nc_s))

    ! Turbulence quantities connected with heat transfer
    if(heat_transfer) then
      call Backup_Mod_Read_Variable(fh,d,vc, 't2',       tur % t2)
      call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'p_t2',     p_t2    (-nb_s:nc_s))
      call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'con_wall', con_wall(-nb_s:nc_s))
    end if
  end if

  !------------------------!
  !   K-eps-zeta-f model   !
  !------------------------!
  if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
     turbulence_model .eq. HYBRID_LES_RANS) then

    ! K, eps, zeta and f22
    call Backup_Mod_Read_Variable(fh,d,vc, 'kin',  tur % kin)
    call Backup_Mod_Read_Variable(fh,d,vc, 'eps',  tur % eps)
    call Backup_Mod_Read_Variable(fh,d,vc, 'zeta', tur % zeta)
    call Backup_Mod_Read_Variable(fh,d,vc, 'f22',  tur % f22)

    ! Other turbulent quantities
    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'p_kin',    p_kin   (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'y_plus',   y_plus  (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'vis_t',    vis_t   (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'vis_wall', vis_wall(-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 't_scale',  tur % t_scale (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'l_scale',  tur % l_scale (-nb_s:nc_s))

    ! Turbulence quantities connected with heat transfer

    if(heat_transfer) then
    call Backup_Mod_Read_Variable(fh,d,vc, 't2',       tur % t2)
    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'p_t2',     p_t2    (-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'con_wall', con_wall(-nb_s:nc_s))
    end if
   
  end if


  !----------------------------!
  !   Reynolds stress models   !
  !----------------------------!
  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then

    ! Reynolds stresses
    call Backup_Mod_Read_Variable(fh,d,vc, 'uu',  tur % uu)
    call Backup_Mod_Read_Variable(fh,d,vc, 'vv',  tur % vv)
    call Backup_Mod_Read_Variable(fh,d,vc, 'ww',  tur % ww)
    call Backup_Mod_Read_Variable(fh,d,vc, 'uv',  tur % uv)
    call Backup_Mod_Read_Variable(fh,d,vc, 'uw',  tur % uw)
    call Backup_Mod_Read_Variable(fh,d,vc, 'vw',  tur % vw)

    ! Epsilon
    call Backup_Mod_Read_Variable(fh,d,vc, 'eps', tur % eps)

    ! F22
    if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
      call Backup_Mod_Read_Variable(fh,d,vc, 'f22',  tur % f22)
    end if

    ! Other turbulent quantities ?
    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'vis_t', vis_t(-nb_s:nc_s))

    ! Turbulence quantities connected with heat transfer
    if(heat_transfer) then
      call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'con_wall', con_wall(-nb_s:nc_s))
    end if
  end if

  !-----------------------------------------!
  !                                         !
  !   Turbulent statistics for all models   !
  !                                         !
  !-----------------------------------------!
  if(turbulence_statistics .and.  &
     time_step > time_step_stat) then

    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'u_mean', tur % u_mean(-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'v_mean', tur % v_mean(-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'w_mean', tur % w_mean(-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'p_mean', tur % p_mean(-nb_s:nc_s))
    if(heat_transfer) then
      call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 't_mean', tur % t_mean(-nb_s:nc_s))
    end if

    ! K and epsilon
    if(turbulence_model .eq. K_EPS) then
      call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'kin_mean',  &
                                          tur % kin_mean(-nb_s:nc_s))
      call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'eps_mean',  &
                                          tur % eps_mean(-nb_s:nc_s))
      if(heat_transfer) then
        call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 't2_mean',  &
                                            tur % t2_mean(-nb_s:nc_s))
      end if
    end if

    ! K-eps-zeta-f and the hybrid model
    if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
       turbulence_model .eq. HYBRID_LES_RANS) then
      call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'kin_mean',  &
                                          tur % kin_mean(-nb_s:nc_s))
      call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'eps_mean',  &
                                          tur % eps_mean(-nb_s:nc_s))
      call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'zeta_mean',  &
                                          tur % zeta_mean(-nb_s:nc_s))
      call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'f22_mean',  &
                                          tur % f22_mean(-nb_s:nc_s))
      if(heat_transfer) then
        call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 't2_mean',  &
                                            tur % t2_mean(-nb_s:nc_s))
      end if
    end if

    ! Reynolds stress models
    if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
       turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
      call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'uu_mean',  &
                                          tur % uu_mean(-nb_s:nc_s))
      call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'vv_mean',  &
                                          tur % vv_mean(-nb_s:nc_s))
      call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'ww_mean',  &
                                          tur % ww_mean(-nb_s:nc_s))
      call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'uv_mean',  &
                                          tur % uv_mean(-nb_s:nc_s))
      call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'uw_mean',  &
                                          tur % uw_mean(-nb_s:nc_s))
      call Backup_Mod_Read_Cell_Bnd(fh, d, vc, 'vw_mean',  &
                                          tur % vw_mean(-nb_s:nc_s))
    end if

    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'uu_res', tur % uu_res(-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'vv_res', tur % vv_res(-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'ww_res', tur % ww_res(-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'uv_res', tur % uv_res(-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'uw_res', tur % uw_res(-nb_s:nc_s))
    call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'vw_res', tur % vw_res(-nb_s:nc_s))

    if(heat_transfer) then
      call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 't2_res', tur % t2_res(-nb_s:nc_s))
      call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'ut_res', tur % ut_res(-nb_s:nc_s))
      call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'vt_res', tur % vt_res(-nb_s:nc_s))
      call Backup_Mod_Read_Cell_Bnd(fh,d,vc, 'wt_res', tur % wt_res(-nb_s:nc_s))
    end if

  end if

  !--------------------------!
  !                          !
  !   Swarm (of particles)   !
  !                          !
  !--------------------------!
  call Backup_Mod_Read_Swarm(fh,d,vc, swr)

  !------------------------------!
  !                              !
  !   User scalars are missing   !
  !                              !
  !------------------------------!

  ! Close backup file
  call Comm_Mod_Close_File(fh)

  end subroutine

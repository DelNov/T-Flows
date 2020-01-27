!==============================================================================!
  subroutine Backup_Mod_Save(fld, swr, tur, mul, &
                             time_step, time_step_stat)
!------------------------------------------------------------------------------!
!   Saves backup files name.backup                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: fld
  type(Swarm_Type),      target :: swr
  type(Turb_Type),       target :: tur
  type(Multiphase_Type), target :: mul
  integer                       :: time_step       ! current time step
  integer                       :: time_step_stat  ! starting step for statist.
!-----------------------------------[Locals]-----------------------------------!
  type(Comm_Type), pointer :: comm
  type(Grid_Type), pointer :: grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: phi
  character(len=80)        :: name_out, name_mean, a_name
  integer                  :: fh, d, vc, sc, ua
!==============================================================================!

  call Cpu_Timer_Mod_Start('Backup_Mode_Save')

  ! Take aliases
  grid => fld % pnt_grid
  bulk => fld % bulk
  comm => grid % comm

  ! Name backup file
  call File_Mod_Set_Name(name_out, time_step=time_step, extension='.backup')

  ! Open backup file
  call Comm_Mod_Open_File_Write(fh, name_out)

  ! Create new types
  call Comm_Mod_Create_New_Types(grid % comm)

  ! Initialize displacement
  d = 0

  ! Intialize number of stored variables
  vc = 0

  !-----------------------------------------------------------------------!
  !   Save cell-centre coordinates.  Could be useful for interpolations   !
  !-----------------------------------------------------------------------!
  call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'x_coords',  &
                                 grid % xc(-comm % nb_s:comm % nc_s))
  call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'y_coords',  &
                                 grid % yc(-comm % nb_s:comm % nc_s))
  call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'z_coords',  &
                                 grid % zc(-comm % nb_s:comm % nc_s))

  !---------------!
  !               !
  !   Save data   !
  !               !
  !---------------!

  ! Time step
  call Backup_Mod_Write_Int(fh, d, vc, 'time_step', time_step)

  ! Number of processors
  call Backup_Mod_Write_Int(fh, d, vc, 'n_proc', n_proc)

  ! Bulk flows and pressure drops in each direction
  call Backup_Mod_Write_Real(fh, d, vc, 'bulk_flux_x',   bulk % flux_x)
  call Backup_Mod_Write_Real(fh, d, vc, 'bulk_flux_y',   bulk % flux_y)
  call Backup_Mod_Write_Real(fh, d, vc, 'bulk_flux_z',   bulk % flux_z)
  call Backup_Mod_Write_Real(fh, d, vc, 'bulk_flux_x_o', bulk % flux_x_o)
  call Backup_Mod_Write_Real(fh, d, vc, 'bulk_flux_y_o', bulk % flux_y_o)
  call Backup_Mod_Write_Real(fh, d, vc, 'bulk_flux_z_o', bulk % flux_z_o)
  call Backup_Mod_Write_Real(fh, d, vc, 'bulk_p_drop_x', bulk % p_drop_x)
  call Backup_Mod_Write_Real(fh, d, vc, 'bulk_p_drop_y', bulk % p_drop_y)
  call Backup_Mod_Write_Real(fh, d, vc, 'bulk_p_drop_z', bulk % p_drop_z)

  !----------------------------!
  !                            !
  !   Navier-Stokes equation   !
  !                            !
  !----------------------------!

  !--------------!
  !   Velocity   !
  !--------------!
  call Backup_Mod_Write_Variable(fh, d, vc, 'u_velocity', fld % u)
  call Backup_Mod_Write_Variable(fh, d, vc, 'v_velocity', fld % v)
  call Backup_Mod_Write_Variable(fh, d, vc, 'w_velocity', fld % w)

  !--------------------------------------!
  !   Pressure and pressure correction   !
  !--------------------------------------!
  call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'press',       &
                                 fld %  p % n(-comm % nb_s:comm % nc_s))
  call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'press_corr',  &
                                 fld % pp % n(-comm % nb_s:comm % nc_s))

  !---------------------!
  !   Mass flow rates   !
  !---------------------!
  call Backup_Mod_Write_Face(grid % comm, fh, d, vc, grid, 'face_flux_00',  &
                             fld % m_flux % n, correct_sign = .true.)

  !--------------!
  !              !
  !   Etnhalpy   !
  !              !
  !--------------!
  if(heat_transfer) then
    call Backup_Mod_Write_Variable(fh, d, vc, 'temp', fld % t)
  end if

  !--------------!
  !              !
  !  Multiphase  !
  !              !
  !--------------!
  if(multiphase_model .eq. VOLUME_OF_FLUID) then
    call Backup_Mod_Write_Variable(fh, d, vc, 'vof', mul % vof)
    call Backup_Mod_Write_Face(grid % comm, fh, d, vc, grid, 'face_dens_00',  &
                               fld % density_f)
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
    call Backup_Mod_Write_Variable(fh, d, vc, 'kin', tur % kin)
    call Backup_Mod_Write_Variable(fh, d, vc, 'eps', tur % eps)

    ! Other turbulent quantities
    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'p_kin',   &
                                   tur % p_kin (-comm % nb_s:comm % nc_s))
    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'y_plus',  &
                                   tur % y_plus(-comm % nb_s:comm % nc_s))
    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'vis_t',   &
                                   tur % vis_t (-comm % nb_s:comm % nc_s))
    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'vis_w',   &
                                   tur % vis_w (-comm % nb_s:comm % nc_s))

    ! Turbulence quantities connected with heat transfer
    if(heat_transfer) then
      call Backup_Mod_Write_Variable(fh, d, vc, 't2',    tur % t2)
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'p_t2',   &
                                     tur % p_t2 (-comm % nb_s:comm % nc_s))
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'con_w',  &
                                     tur % con_w(-comm % nb_s:comm % nc_s))
    end if

  end if

  !------------------------!
  !   K-eps-zeta-f model   !
  !------------------------!
  if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
     turbulence_model .eq. HYBRID_LES_RANS) then

    ! K, eps, zeta and f22
    call Backup_Mod_Write_Variable(fh, d, vc, 'kin',  tur % kin)
    call Backup_Mod_Write_Variable(fh, d, vc, 'eps',  tur % eps)
    call Backup_Mod_Write_Variable(fh, d, vc, 'zeta', tur % zeta)
    call Backup_Mod_Write_Variable(fh, d, vc, 'f22',  tur % f22)

    ! Other turbulent quantities
    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc,'p_kin',    &
                                   tur % p_kin  (-comm % nb_s:comm % nc_s))
    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc,'y_plus',   &
                                   tur % y_plus (-comm % nb_s:comm % nc_s))
    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc,'vis_t',    &
                                   tur % vis_t  (-comm % nb_s:comm % nc_s))
    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc,'vis_w',    &
                                   tur % vis_w  (-comm % nb_s:comm % nc_s))
    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc,'t_scale',  &
                                   tur % t_scale(-comm % nb_s:comm % nc_s))
    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc,'l_scale',  &
                                   tur % l_scale(-comm % nb_s:comm % nc_s))

    if(heat_transfer) then
      call Backup_Mod_Write_Variable(fh, d, vc, 't2',    tur % t2)
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'p_t2',   &
                                     tur % p_t2 (-comm % nb_s:comm % nc_s))
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'con_w',  &
                                     tur % con_w(-comm % nb_s:comm % nc_s))
    end if

  end if

  !----------------------------!
  !   Reynolds stress models   !
  !----------------------------!
  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then

    ! Reynolds stresses
    call Backup_Mod_Write_Variable(fh, d, vc, 'uu',  tur % uu)
    call Backup_Mod_Write_Variable(fh, d, vc, 'vv',  tur % vv)
    call Backup_Mod_Write_Variable(fh, d, vc, 'ww',  tur % ww)
    call Backup_Mod_Write_Variable(fh, d, vc, 'uv',  tur % uv)
    call Backup_Mod_Write_Variable(fh, d, vc, 'uw',  tur % uw)
    call Backup_Mod_Write_Variable(fh, d, vc, 'vw',  tur % vw)

    ! Epsilon
    call Backup_Mod_Write_Variable(fh, d, vc, 'eps', tur % eps)

    ! F22
    if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
      call Backup_Mod_Write_Variable(fh, d, vc, 'f22',  tur % f22)
    end if

    ! Other turbulent quantities 
    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'vis_t',  &
                                   tur % vis_t(-comm % nb_s:comm % nc_s))
    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'vis_w',  &
                                   tur % vis_w(-comm % nb_s:comm % nc_s))

    ! Turbulence quantities connected with heat transfer
    if(heat_transfer) then
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc,'con_w', &
                                     tur % con_w(-comm % nb_s:comm % nc_s))
    end if
  end if

  !------------------!
  !   Save scalars   !
  !------------------!
  do sc = 1, fld % n_scalars
    phi => fld % scalar(sc)
    call Backup_Mod_Write_Variable(fh, d, vc, phi % name, phi)
  end do

  !-----------------------------------------!
  !                                         !
  !   Turbulent statistics for all models   !
  !                                         !
  !-----------------------------------------!
  if(turbulence_statistics .and.  &
     time_step > time_step_stat) then

    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'u_mean',  &
                                   tur % u_mean(-comm % nb_s:comm % nc_s))
    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'v_mean',  &
                                   tur % v_mean(-comm % nb_s:comm % nc_s))
    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'w_mean',  &
                                   tur % w_mean(-comm % nb_s:comm % nc_s))
    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'p_mean',  &
                                   tur % p_mean(-comm % nb_s:comm % nc_s))
    if(heat_transfer) then
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 't_mean', &
                                     tur % t_mean(-comm % nb_s:comm % nc_s))
    end if

    ! K and epsilon
    if(turbulence_model .eq. K_EPS) then
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'kin_mean',  &
                                     tur % kin_mean(-comm % nb_s:comm % nc_s))
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'eps_mean',  &
                                     tur % eps_mean(-comm % nb_s:comm % nc_s))
      if(heat_transfer) then
        call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 't2_mean',  &
                                       tur % t2_mean(-comm % nb_s:comm % nc_s))
      end if
    end if

    ! K-eps-zeta-f and the hybrid model
    if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
       turbulence_model .eq. HYBRID_LES_RANS) then
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'kin_mean',  &
                                     tur % kin_mean (-comm % nb_s:comm % nc_s))
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'eps_mean',  &
                                     tur % eps_mean (-comm % nb_s:comm % nc_s))
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'zeta_mean',  &
                                     tur % zeta_mean(-comm % nb_s:comm % nc_s))
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'f22_mean',  &
                                     tur % f22_mean (-comm % nb_s:comm % nc_s))
      if(heat_transfer) then
        call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 't2_mean',  &
                                       tur % t2_mean(-comm % nb_s:comm % nc_s))
      end if
    end if

    ! Reynolds stress models
    if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
       turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'uu_mean',  &
                                     tur % uu_mean(-comm % nb_s:comm % nc_s))
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'vv_mean',  &
                                     tur % vv_mean(-comm % nb_s:comm % nc_s))
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'ww_mean',  &
                                     tur % ww_mean(-comm % nb_s:comm % nc_s))
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'uv_mean',  &
                                     tur % uv_mean(-comm % nb_s:comm % nc_s))
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'uw_mean',  &
                                     tur % uw_mean(-comm % nb_s:comm % nc_s))
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'vw_mean',  &
                                     tur % vw_mean(-comm % nb_s:comm % nc_s))
    end if

    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'uu_res',  &
                                   tur % uu_res(-comm % nb_s:comm % nc_s))
    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'vv_res',  &
                                   tur % vv_res(-comm % nb_s:comm % nc_s))
    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'ww_res',  &
                                   tur % ww_res(-comm % nb_s:comm % nc_s))
    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'uv_res',  &
                                   tur % uv_res(-comm % nb_s:comm % nc_s))
    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'uw_res',  &
                                   tur % uw_res(-comm % nb_s:comm % nc_s))
    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'vw_res',  &
                                   tur % vw_res(-comm % nb_s:comm % nc_s))

    if(heat_transfer) then
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 't2_res',  &
                                     tur % t2_res(-comm % nb_s:comm % nc_s))
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'ut_res',  &
                                     tur % ut_res(-comm % nb_s:comm % nc_s))
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'vt_res',  &
                                     tur % vt_res(-comm % nb_s:comm % nc_s))
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, 'wt_res',  &
                                     tur % wt_res(-comm % nb_s:comm % nc_s))
    end if

    ! Scalars
    do sc = 1, fld % n_scalars
      phi => fld % scalar(sc)
      name_mean = phi % name
      name_mean(5:9) = '_mean'
      call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, name_mean,  &
                            tur % scalar_mean(sc, -comm % nb_s:comm % nc_s))
    end do

  end if

  !--------------------------!
  !                          !
  !   Swarm (of particles)   !
  !                          !
  !--------------------------!
  call Backup_Mod_Write_Swarm(fh, d, vc, swr)

  !-----------------!
  !                 !
  !   User arrays   !
  !                 !
  !-----------------!

  do ua = 1, n_user_arrays
    a_name = 'A_??'
    write(a_name(3:4),'(I2.2)') ua
    call Backup_Mod_Write_Cell_Bnd(comm, fh, d, vc, a_name,  &
                                   user_array(ua,-comm % nb_s:comm % nc_s))
  end do

  ! Variable count (store +1 to count its own self)
  call Backup_Mod_Write_Int(fh, d, vc, 'variable_count', vc + 1)

  if(this_proc < 2) then
    print *, '# Wrote ', vc, ' variables!'
  end if

  ! Close backup file
  call Comm_Mod_Close_File(fh)

  call Cpu_Timer_Mod_Stop('Backup_Mode_Save')

  end subroutine

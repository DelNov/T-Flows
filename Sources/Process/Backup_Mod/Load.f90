!==============================================================================!
  subroutine Backup_Mod_Load(Fld, Swr, tur, Vof, time, time_step, backup)
!------------------------------------------------------------------------------!
!   Loads backup files name.backup                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Fld
  type(Swarm_Type), target :: Swr
  type(Turb_Type),  target :: tur
  type(Vof_Type),   target :: Vof
  real                     :: time            ! time of simulation
  integer                  :: time_step       ! current time step
  logical                  :: backup, present
!-----------------------------------[Locals]-----------------------------------!
  type(Comm_Type), pointer :: Comm
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: phi
  character(SL)            :: name_in, answer, name_mean, a_name
  integer                  :: fh, d, vc, sc, ua
!==============================================================================!

  call Cpu_Timer % Start('Backup_Mod_Load')

  ! Take aliases
  Grid => Fld % pnt_grid
  bulk => Fld % bulk
  Comm => Grid % Comm

  ! Full name is specified in control file
  call Control_Mod_Load_Backup_Name(name_in)

  answer = name_in
  call To_Upper_Case(answer)

  backup = .true.
  if(answer .eq. 'SKIP') then
    backup = .false.
    call Cpu_Timer % Stop('Backup_Mod_Load')
    return
  end if

  inquire(file=trim(name_in), exist=present )
  if(.not.present) then
    if(this_proc < 2) then
      print *, "# ERROR!  Backup file ", trim(name_in), " was not found."
      print *, "# Exiting!"
    end if
    call Comm_Mod_End
    stop
  end if

  ! Open backup file
  call Comm % Open_File_Read(fh, name_in)

  ! Create new types
  call Comm % Create_New_Types()

  ! Initialize displacement and variable count
  d  = 0
  vc = 0

  !---------------------------------------------!
  !   Variable count - important for checking   !
  !---------------------------------------------!
  call Backup_Mod_Read_Int(Comm, fh, d, 2048, 'variable_count', vc)
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
  ! call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'x_coords', Grid % xc)
  ! call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'y_coords', Grid % yc)
  ! call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'z_coords', Grid % zc)

  ! Time step
  call Backup_Mod_Read_Int(Comm, fh, d, vc, 'time_step', time_step)

  ! Simulation time
  call Backup_Mod_Read_Real(Comm, fh, d, vc, 'time', time)

  ! Bulk flows and pressure drops in each direction
  call Backup_Mod_Read_Real(Comm, fh, d, vc, 'bulk_flux_x',   bulk % flux_x)
  call Backup_Mod_Read_Real(Comm, fh, d, vc, 'bulk_flux_y',   bulk % flux_y)
  call Backup_Mod_Read_Real(Comm, fh, d, vc, 'bulk_flux_z',   bulk % flux_z)
  call Backup_Mod_Read_Real(Comm, fh, d, vc, 'bulk_flux_x_o', bulk % flux_x_o)
  call Backup_Mod_Read_Real(Comm, fh, d, vc, 'bulk_flux_y_o', bulk % flux_y_o)
  call Backup_Mod_Read_Real(Comm, fh, d, vc, 'bulk_flux_z_o', bulk % flux_z_o)
  call Backup_Mod_Read_Real(Comm, fh, d, vc, 'bulk_p_drop_x', bulk % p_drop_x)
  call Backup_Mod_Read_Real(Comm, fh, d, vc, 'bulk_p_drop_y', bulk % p_drop_y)
  call Backup_Mod_Read_Real(Comm, fh, d, vc, 'bulk_p_drop_z', bulk % p_drop_z)

  !----------------------------!
  !                            !
  !   Navier-Stokes equation   !
  !                            !
  !----------------------------!

  !--------------!
  !   Velocity   !
  !--------------!
  call Backup_Mod_Read_Variable(fh, d, vc, 'u_velocity', Fld, Fld % u)
  call Backup_Mod_Read_Variable(fh, d, vc, 'v_velocity', Fld, Fld % v)
  call Backup_Mod_Read_Variable(fh, d, vc, 'w_velocity', Fld, Fld % w)

  !------------------------------------------------------!
  !   Pressure, its gradients, and pressure correction   !
  !------------------------------------------------------!
  call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'press',      Fld % p % n)
  call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'press_x',    Fld % p % x)
  call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'press_y',    Fld % p % y)
  call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'press_z',    Fld % p % z)
  call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'press_corr', Fld % pp % n)
  call Fld % Grad_Pressure(Fld % pp)

  !-------------------!
  !   Volume fluxes   ! -> don't use for the time being, too much trouble
  !-------------------!

  !----------------------------------!
  !   Does it have pressure outlet   !
  !----------------------------------!
  ! Update on July 17, 2021: I have some reservations about this part, since
  ! there was another bug fix when computing fluxes in the meanwhile (check:
  ! 90f77a1c8bd4ca05330a4435ed6321782ef00199).  This balancing also caused a
  ! bug when loading backup file (also check "Compute_Pressure" as well as
  ! "Initialize_Variables" procedures)
  call Backup_Mod_Read_Log(Comm, fh, d, vc, 'has_pressure_outlet',  &
                           Fld % has_pressure_outlet)

  !--------------!
  !              !
  !   Enthalpy   !
  !              !
  !--------------!
  if(Fld % heat_transfer) then
    call Backup_Mod_Read_Variable(fh, d, vc, 'temp', Fld, Fld % t)
  end if

  !--------------!
  !              !
  !   Roughness  !
  !              !
  !--------------!
  if(tur % rough_walls) then
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc,'id_zone',  &
                                   tur % id_zone)
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc,'z_o_f',    &
                                   tur % z_o_f)
  end if

  !--------------!
  !              !
  !  Multiphase  !
  !              !
  !--------------!
  if(Fld % with_interface) then
    call Backup_Mod_Read_Variable(fh, d, vc, 'vof_fun', Fld, Vof % fun)
  end if

  !-----------------------!
  !                       !
  !   Turbulence models   !
  !                       !
  !-----------------------!

  !-----------------!
  !   K-eps model   !
  !-----------------!
  if(tur % model .eq. K_EPS) then

    ! K and epsilon
    call Backup_Mod_Read_Variable(fh, d, vc, 'kin', Fld, tur % kin)
    call Backup_Mod_Read_Variable(fh, d, vc, 'eps', Fld, tur % eps)

    ! Other turbulent quantities
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'p_kin',  tur % p_kin )
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'y_plus', tur % y_plus)
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'vis_t',  tur % vis_t )
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'vis_w',  tur % vis_w )

    ! Turbulence quantities connected with heat transfer
    if(Fld % heat_transfer) then
      call Backup_Mod_Read_Variable(fh, d, vc, 't2', Fld, tur % t2)
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'p_t2',  tur % p_t2 )
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'con_w', tur % con_w)
    end if
  end if

  !------------------------!
  !   K-eps-zeta-f model   !
  !------------------------!
  if(tur % model .eq. K_EPS_ZETA_F .or.  &
     tur % model .eq. HYBRID_LES_RANS) then

    ! K, eps, zeta and f22
    call Backup_Mod_Read_Variable(fh, d, vc, 'kin',  Fld, tur % kin)
    call Backup_Mod_Read_Variable(fh, d, vc, 'eps',  Fld, tur % eps)
    call Backup_Mod_Read_Variable(fh, d, vc, 'zeta', Fld, tur % zeta)
    call Backup_Mod_Read_Variable(fh, d, vc, 'f22',  Fld, tur % f22)

    ! Other turbulent quantities
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc,'p_kin',   tur % p_kin  )
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc,'y_plus',  tur % y_plus )
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc,'vis_t',   tur % vis_t  )
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc,'vis_w',   tur % vis_w  )
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc,'t_scale', tur % t_scale)
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc,'l_scale', tur % l_scale)

    ! Turbulence quantities connected with heat transfer

    if(Fld % heat_transfer) then
      call Backup_Mod_Read_Variable(fh, d, vc, 't2', Fld, tur % t2)
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'p_t2',  tur % p_t2 )
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'con_w', tur % con_w)
    end if

  end if

  !----------------------------!
  !   Reynolds stress models   !
  !----------------------------!
  if(tur % model .eq. RSM_MANCEAU_HANJALIC .or.  &
     tur % model .eq. RSM_HANJALIC_JAKIRLIC) then

    ! Reynolds stresses
    call Backup_Mod_Read_Variable(fh, d, vc, 'uu', Fld, tur % uu)
    call Backup_Mod_Read_Variable(fh, d, vc, 'vv', Fld, tur % vv)
    call Backup_Mod_Read_Variable(fh, d, vc, 'ww', Fld, tur % ww)
    call Backup_Mod_Read_Variable(fh, d, vc, 'uv', Fld, tur % uv)
    call Backup_Mod_Read_Variable(fh, d, vc, 'uw', Fld, tur % uw)
    call Backup_Mod_Read_Variable(fh, d, vc, 'vw', Fld, tur % vw)

    ! Epsilon
    call Backup_Mod_Read_Variable(fh, d, vc, 'eps', Fld, tur % eps)

    ! F22
    if(tur % model .eq. RSM_MANCEAU_HANJALIC) then
      call Backup_Mod_Read_Variable(fh, d, vc, 'f22', Fld, tur % f22)
    end if

    ! Other turbulent quantities ?
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'vis_t', tur % vis_t)
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'vis_w', tur % vis_w)

    ! Turbulence quantities connected with heat transfer
    if(Fld % heat_transfer) then
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'con_w', tur % con_w)
    end if
  end if

  !------------------!
  !   Load scalars   !
  !------------------!
  do sc = 1, Fld % n_scalars
    phi => Fld % scalar(sc)
    call Backup_Mod_Read_Variable(fh, d, vc, phi % name, Fld, phi)
  end do

  !-----------------------------------------!
  !                                         !
  !   Turbulent statistics for all models   !
  !                                         !
  !-----------------------------------------!
  if(tur % statistics) then

    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'u_mean', tur % u_mean)
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'v_mean', tur % v_mean)
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'w_mean', tur % w_mean)
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'p_mean', tur % p_mean)
    if(Fld % heat_transfer) then
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 't_mean', tur % t_mean)
    end if

    ! K and epsilon
    if(tur % model .eq. K_EPS) then
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'kin_mean',  &
                                     tur % kin_mean)
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'eps_mean',  &
                                     tur % eps_mean)
      if(Fld % heat_transfer) then
        call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 't2_mean',  &
                                       tur % t2_mean)
        call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'ut_mean',  &
                                       tur % ut_mean)
        call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'vt_mean',  &
                                       tur % vt_mean)
        call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'wt_mean',  &
                                       tur % wt_mean)
      end if
    end if

    ! K-eps-zeta-f and the hybrid model
    if(tur % model .eq. K_EPS_ZETA_F .or.  &
       tur % model .eq. HYBRID_LES_RANS) then
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'kin_mean',  &
                                     tur % kin_mean)
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'eps_mean',  &
                                     tur % eps_mean)
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'zeta_mean',  &
                                     tur % zeta_mean)
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'f22_mean',  &
                                     tur % f22_mean)
      if(Fld % heat_transfer) then
        call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 't2_mean',  &
                                       tur % t2_mean)
        call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'ut_mean',  &
                                       tur % ut_mean)
        call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'vt_mean',  &
                                       tur % vt_mean)
        call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'wt_mean',  &
                                       tur % wt_mean)
      end if
    end if

    ! Reynolds stress models
    if(tur % model .eq. RSM_MANCEAU_HANJALIC .or.  &
       tur % model .eq. RSM_HANJALIC_JAKIRLIC) then
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'uu_mean',  &
                                     tur % uu_mean)
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'vv_mean',  &
                                     tur % vv_mean)
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'ww_mean',  &
                                     tur % ww_mean)
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'uv_mean',  &
                                     tur % uv_mean)
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'uw_mean',  &
                                     tur % uw_mean)
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'vw_mean',  &
                                     tur % vw_mean)
    end if

    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'uu_res', tur % uu_res)
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'vv_res', tur % vv_res)
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'ww_res', tur % ww_res)
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'uv_res', tur % uv_res)
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'uw_res', tur % uw_res)
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'vw_res', tur % vw_res)

    if(Fld % heat_transfer) then
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 't2_res', tur % t2_res)
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'ut_res', tur % ut_res)
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'vt_res', tur % vt_res)
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'wt_res', tur % wt_res)
    end if

    ! Scalars
    do sc = 1, Fld % n_scalars
      phi => Fld % scalar(sc)
      name_mean = phi % name
      name_mean(5:9) = '_mean'
      call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, name_mean,  &
                                     tur % scalar_mean(sc, :))
    end do

  end if

  !--------------------------!
  !                          !
  !   Swarm (of particles)   !
  !                          !
  !--------------------------!
  if(Fld % with_particles) then
    call Backup_Mod_Read_Swarm(fh, d, vc, Swr)
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'n_deposited',      &
                                   Swr % n_deposited)
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, 'n_reflected',      &
                                   Swr % n_reflected)
  end if

  !-----------------!
  !                 !
  !   User arrays   !
  !                 !
  !-----------------!

  do ua = 1, Grid % n_user_arrays
    a_name = 'A_00'
    write(a_name(3:4),'(I2.2)') ua
    call Backup_Mod_Read_Cell_Real(Grid, fh, d, vc, a_name,  &
                                   Grid % user_array(ua, :))
  end do

  ! Close backup file
  call Comm % Close_File(fh)

  call Cpu_Timer % Stop('Backup_Mod_Load')

  end subroutine

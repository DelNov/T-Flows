!==============================================================================!
  subroutine Load(Backup, Flow, Turb, Vof, Swarm, bckp)
!------------------------------------------------------------------------------!
!   Loads backup files name.backup                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Backup_Type)       :: Backup
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  logical                  :: bckp, present
!-----------------------------------[Locals]-----------------------------------!
  type(Comm_Type), pointer :: Comm
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: phi
  character(SL)            :: name_in, answer, name_mean
  integer                  :: vc, sc, ts
  integer(DP)              :: d
  real                     :: st  ! saved time, simulation time
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  bulk => Flow % bulk
  Comm => Grid % Comm

  ! Full name is specified in control file
  call Control % Load_Backup_Name(name_in)

  answer = name_in
  call String % To_Upper_Case(answer)

  bckp = .true.
  if(answer .eq. 'SKIP') then
    bckp = .false.
    return
  end if

  ! OK, you are rading a backup, you may as well time it
  call Profiler % Start('Backup_Mod_Load')

  inquire(file=trim(name_in), exist=present )
  if(.not.present) then
    call Message % Error(80,                                                &
               'Backup file '//trim(name_in)//' was not found.  Exiting!',  &
               file=__FILE__, line=__LINE__, one_proc=.true.)
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
  call Backup % Load_Int(Comm, d, 2048, 'variable_count', vc)
  if(vc .eq. 0) vc = 2048  ! for backward compatibility

  if(First_Proc()) then
    print *, "# Backup file holds ", vc, " variables."
  end if

  !---------------!
  !               !
  !   Load data   !
  !               !
  !---------------!

  !---------------------------------!
  !   Related to time integration   !
  !---------------------------------!

  ! Time step
  call Backup % Load_Int(Comm, d, vc, 'time_step', ts)
  call Time % Set_First_Dt(ts)

  ! Simulation time
  call Backup % Load_Real(Comm, d, vc, 'time', st)
  call Time % Set_Time(st)

  !-----------------------------------------------------!
  !   Bulk flows and pressure drops in each direction   !
  !-----------------------------------------------------!
  call Backup % Load_Real(Comm, d, vc, 'bulk_flux_x',   bulk % flux_x)
  call Backup % Load_Real(Comm, d, vc, 'bulk_flux_y',   bulk % flux_y)
  call Backup % Load_Real(Comm, d, vc, 'bulk_flux_z',   bulk % flux_z)
  call Backup % Load_Real(Comm, d, vc, 'bulk_flux_x_o', bulk % flux_x_o)
  call Backup % Load_Real(Comm, d, vc, 'bulk_flux_y_o', bulk % flux_y_o)
  call Backup % Load_Real(Comm, d, vc, 'bulk_flux_z_o', bulk % flux_z_o)
  call Backup % Load_Real(Comm, d, vc, 'bulk_p_drop_x', bulk % p_drop_x)
  call Backup % Load_Real(Comm, d, vc, 'bulk_p_drop_y', bulk % p_drop_y)
  call Backup % Load_Real(Comm, d, vc, 'bulk_p_drop_z', bulk % p_drop_z)

  !----------------------------!
  !                            !
  !   Navier-Stokes equation   !
  !                            !
  !----------------------------!

  !--------------!
  !   Velocity   !
  !--------------!
  call Backup % Load_Variable(d, vc, 'u_velocity', Flow, Flow % u)
  call Backup % Load_Variable(d, vc, 'v_velocity', Flow, Flow % v)
  call Backup % Load_Variable(d, vc, 'w_velocity', Flow, Flow % w)

  !------------------------------------------------------!
  !   Pressure, its gradients, and pressure correction   !
  !------------------------------------------------------!
  call Backup % Load_Cell_Real(Grid, d, vc, 'press',      Flow % p % n)
  call Backup % Load_Cell_Real(Grid, d, vc, 'press_x',    Flow % p % x)
  call Backup % Load_Cell_Real(Grid, d, vc, 'press_y',    Flow % p % y)
  call Backup % Load_Cell_Real(Grid, d, vc, 'press_z',    Flow % p % z)
  call Backup % Load_Cell_Real(Grid, d, vc, 'press_corr', Flow % pp % n)
  call Flow % Grad_Pressure(Flow % pp)

  !-------------------!
  !   Volume fluxes   ! -> don't use for the time being, too much trouble
  !-------------------!

  !-------------------------!
  !   Does it have outlet   !
  !-------------------------!
  ! Update on July 17, 2021: I have some reservations about this part, since
  ! there was another bug fix when computing fluxes in the meanwhile (check:
  ! 90f77a1c8bd4ca05330a4435ed6321782ef00199).  This balancing also caused a
  ! bug when loading backup file (also check "Compute_Pressure" as well as
  ! "Initialize_Variables" procedures)
  !
  ! Update on February 27, 2022: I have also added "has_outflow_boundary"
  ! to be able to tell PETSc if matrix for pressure is singular
  !
  ! Update on June 2, 2022: Unified all outlet boundaries into one
  ! to be able to tell PETSc if matrix for pressure is singular
  call Backup % Load_Log(Comm, d, vc, 'has_pressure', Flow % has_pressure)

  !--------------!
  !              !
  !   Enthalpy   !
  !              !
  !--------------!
  if(Flow % heat_transfer) then
    call Backup % Load_Variable(d, vc, 'temp', Flow, Flow % t)
  end if

  !--------------!
  !              !
  !  Multiphase  !
  !              !
  !--------------!
  if(Flow % with_interface) then
    call Backup % Load_Variable(d, vc, 'vof_fun', Flow, Vof % fun)
  end if

  !-----------------------!
  !                       !
  !   Turbulence models   !
  !                       !
  !-----------------------!

  !-----------------!
  !   K-eps model   !
  !-----------------!
  if(Turb % model .eq. K_EPS) then

    ! K and epsilon
    call Backup % Load_Variable(d, vc, 'kin', Flow, Turb % kin)
    call Backup % Load_Variable(d, vc, 'eps', Flow, Turb % eps)

    ! Other turbulent quantities
    call Backup % Load_Cell_Real(Grid, d, vc, 'p_kin',  Turb % p_kin )
    call Backup % Load_Cell_Real(Grid, d, vc, 'y_plus', Turb % y_plus)
    call Backup % Load_Cell_Real(Grid, d, vc, 'vis_t',  Turb % vis_t )
    call Backup % Load_Cell_Real(Grid, d, vc, 'vis_w',  Turb % vis_w )

    ! Turbulence quantities connected with heat transfer
    if(Flow % heat_transfer) then
      call Backup % Load_Variable(d, vc, 't2', Flow, Turb % t2)
      call Backup % Load_Cell_Real(Grid, d, vc, 'p_t2',  Turb % p_t2 )
      call Backup % Load_Cell_Real(Grid, d, vc, 'con_w', Turb % con_w)
    end if
  end if

  !------------------------!
  !   K-eps-zeta-f model   !
  !------------------------!
  if(Turb % model .eq. K_EPS_ZETA_F .or.  &
     Turb % model .eq. HYBRID_LES_RANS) then

    ! K, eps, zeta and f22
    call Backup % Load_Variable(d, vc, 'kin',  Flow, Turb % kin)
    call Backup % Load_Variable(d, vc, 'eps',  Flow, Turb % eps)
    call Backup % Load_Variable(d, vc, 'zeta', Flow, Turb % zeta)
    call Backup % Load_Variable(d, vc, 'f22',  Flow, Turb % f22)

    ! Other turbulent quantities
    call Backup % Load_Cell_Real(Grid, d, vc,'p_kin',   Turb % p_kin  )
    call Backup % Load_Cell_Real(Grid, d, vc,'y_plus',  Turb % y_plus )
    call Backup % Load_Cell_Real(Grid, d, vc,'vis_t',   Turb % vis_t  )
    call Backup % Load_Cell_Real(Grid, d, vc,'vis_w',   Turb % vis_w  )
    call Backup % Load_Cell_Real(Grid, d, vc,'t_scale', Turb % t_scale)
    call Backup % Load_Cell_Real(Grid, d, vc,'l_scale', Turb % l_scale)

    ! Turbulence quantities connected with heat transfer

    if(Flow % heat_transfer) then
      call Backup % Load_Variable(d, vc, 't2', Flow, Turb % t2)
      call Backup % Load_Cell_Real(Grid, d, vc, 'p_t2',  Turb % p_t2 )
      call Backup % Load_Cell_Real(Grid, d, vc, 'con_w', Turb % con_w)
    end if

  end if

  !-------------------------!
  !   K-omega-sst model     !
  !-------------------------!
  if(Turb % model .eq. K_OMEGA_SST) then
    call Backup % Load_Variable(d, vc, 'kin', Flow, Turb % kin)
    call Backup % Load_Variable(d, vc, 'omg', Flow, Turb % omega)

    call Backup % Load_Cell_Real(Grid, d, vc,'p_kin',   Turb % p_kin )
    call Backup % Load_Cell_Real(Grid, d, vc,'y_plus',  Turb % y_plus)
    call Backup % Load_Cell_Real(Grid, d, vc,'vis_t',   Turb % vis_t )
    call Backup % Load_Cell_Real(Grid, d, vc,'vis_w',   Turb % vis_w )
    call Backup % Load_Cell_Real(Grid, d, vc,'t_scale', Turb % t_scale)
  end if

  !-----------------!
  !   S-A model     !
  !-----------------!
  if(Turb % model .eq. SPALART_ALLMARAS) then

    ! K and epsilon
    call Backup % Load_Variable(d, vc, 'vis', Flow, Turb % vis)

    ! Other turbulent quantities
    call Backup % Load_Cell_Real(Grid, d, vc, 'y_plus', Turb % y_plus)
    call Backup % Load_Cell_Real(Grid, d, vc, 'vis_t',  Turb % vis_t )
    call Backup % Load_Cell_Real(Grid, d, vc, 'vis_w',  Turb % vis_w )

    ! Turbulence quantities connected with heat transfer
    if(Flow % heat_transfer) then
      call Backup % Load_Cell_Real(Grid, d, vc, 'con_w', Turb % con_w)
    end if

  end if

  !--------------!
  !   Roughness  !
  !--------------!
  call Backup % Load_Cell_Real(Grid, d, vc, 'z_o', Turb % z_o)

  !------------------!
  !   Load scalars   !
  !------------------!
  do sc = 1, Flow % n_scalars
    phi => Flow % scalar(sc)
    call Backup % Load_Variable(d, vc, phi % name, Flow, phi)
  end do

  !-----------------------------------------!
  !                                         !
  !   Turbulent statistics for all models   !
  !                                         !
  !-----------------------------------------!
  if(Turb % statistics) then

    call Backup % Load_Cell_Real(Grid, d, vc, 'u_mean', Turb % u_mean)
    call Backup % Load_Cell_Real(Grid, d, vc, 'v_mean', Turb % v_mean)
    call Backup % Load_Cell_Real(Grid, d, vc, 'w_mean', Turb % w_mean)
    call Backup % Load_Cell_Real(Grid, d, vc, 'p_mean', Turb % p_mean)
    if(Flow % heat_transfer) then
      call Backup % Load_Cell_Real(Grid, d, vc, 't_mean', Turb % t_mean)
      call Backup % Load_Cell_Real(Grid, d, vc, 'q_mean', Turb % q_mean)
    end if

    ! K and epsilon
    if(Turb % model .eq. K_EPS) then
      call Backup % Load_Cell_Real(Grid, d, vc, 'kin_mean', Turb % kin_mean)
      call Backup % Load_Cell_Real(Grid, d, vc, 'eps_mean', Turb % eps_mean)
      if(Flow % heat_transfer) then
        call Backup % Load_Cell_Real(Grid, d, vc, 't2_mean', Turb % t2_mean)
        call Backup % Load_Cell_Real(Grid, d, vc, 'ut_mean', Turb % ut_mean)
        call Backup % Load_Cell_Real(Grid, d, vc, 'vt_mean', Turb % vt_mean)
        call Backup % Load_Cell_Real(Grid, d, vc, 'wt_mean', Turb % wt_mean)
      end if
    end if

    ! K-eps-zeta-f and the hybrid model
    if(Turb % model .eq. K_EPS_ZETA_F .or.  &
       Turb % model .eq. HYBRID_LES_RANS) then
      call Backup % Load_Cell_Real(Grid, d, vc, 'kin_mean',  Turb % kin_mean)
      call Backup % Load_Cell_Real(Grid, d, vc, 'eps_mean',  Turb % eps_mean)
      call Backup % Load_Cell_Real(Grid, d, vc, 'zeta_mean', Turb % zeta_mean)
      call Backup % Load_Cell_Real(Grid, d, vc, 'f22_mean',  Turb % f22_mean)
      if(Flow % heat_transfer) then
        call Backup % Load_Cell_Real(Grid, d, vc, 't2_mean', Turb % t2_mean)
        call Backup % Load_Cell_Real(Grid, d, vc, 'ut_mean', Turb % ut_mean)
        call Backup % Load_Cell_Real(Grid, d, vc, 'vt_mean', Turb % vt_mean)
        call Backup % Load_Cell_Real(Grid, d, vc, 'wt_mean', Turb % wt_mean)
      end if
    end if

    call Backup % Load_Cell_Real(Grid, d, vc, 'uu_res', Turb % uu_res)
    call Backup % Load_Cell_Real(Grid, d, vc, 'vv_res', Turb % vv_res)
    call Backup % Load_Cell_Real(Grid, d, vc, 'ww_res', Turb % ww_res)
    call Backup % Load_Cell_Real(Grid, d, vc, 'uv_res', Turb % uv_res)
    call Backup % Load_Cell_Real(Grid, d, vc, 'uw_res', Turb % uw_res)
    call Backup % Load_Cell_Real(Grid, d, vc, 'vw_res', Turb % vw_res)

    if(Flow % heat_transfer) then
      call Backup % Load_Cell_Real(Grid, d, vc, 't2_res', Turb % t2_res)
      call Backup % Load_Cell_Real(Grid, d, vc, 'ut_res', Turb % ut_res)
      call Backup % Load_Cell_Real(Grid, d, vc, 'vt_res', Turb % vt_res)
      call Backup % Load_Cell_Real(Grid, d, vc, 'wt_res', Turb % wt_res)
    end if

    ! Scalars
    do sc = 1, Flow % n_scalars
      phi => Flow % scalar(sc)
      name_mean = phi % name
      name_mean(5:9) = '_mean'
      call Backup % Load_Cell_Real(Grid, d, vc, name_mean,  &
                                   Turb % scalar_mean(sc, :))
    end do

  end if

  !--------------------------!
  !                          !
  !   Swarm (of particles)   !
  !                          !
  !--------------------------!
  if(Flow % with_particles) then
    call Backup % Load_Swarm(Swarm, d, vc)
    call Backup % Load_Cell_Real(Grid, d,vc,'n_deposited', Swarm % n_deposited)
    call Backup % Load_Cell_Real(Grid, d,vc,'n_reflected', Swarm % n_reflected)
  end if

  ! Close backup file
  call Comm % Close_File(fh)

  call Profiler % Stop('Backup_Mod_Load')

  end subroutine

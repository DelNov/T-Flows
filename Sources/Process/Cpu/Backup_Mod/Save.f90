!==============================================================================!
  subroutine Save(Backup, Flow, Turb, Vof, Swarm, dom)
!------------------------------------------------------------------------------!
!   Saves backup files name.backup                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Backup_Type)        :: Backup
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: Turb
  type(Vof_Type),    target :: Vof
  type(Swarm_Type),  target :: Swarm
  integer,         optional :: dom
!-----------------------------------[Locals]-----------------------------------!
  type(Comm_Type), pointer :: Comm
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: phi
  character(SL)            :: name_out, name_mean
  integer                  :: vc, sc
  integer(DP)              :: d
!==============================================================================!

  call Profiler % Start('Backup_Mod_Save')

  ! Take aliases
  Grid => Flow % pnt_grid
  bulk => Flow % bulk
  Comm => Grid % Comm

  ! Name backup file
  call File % Set_Name(name_out, time_step = Time % Curr_Dt(),  &
                       extension='.backup', domain=dom)

  ! Open backup file
  call Comm % Open_File_Write(fh, name_out)

  ! Create new types
  call Comm % Create_New_Types()

  ! Initialize displacement
  d = 0

  ! Intialize number of stored variables
  vc = 0

  !---------------!
  !               !
  !   Save data   !
  !               !
  !---------------!

  !---------------------------------!
  !   Related to time integration   !
  !---------------------------------!

  ! Time step
  call Backup % Save_Int(Comm, d, vc, 'time_step', Time % Curr_Dt())

  ! Simulation time
  call Backup % Save_Real(Comm, d, vc, 'time', Time % Get_Time())

  ! Number of processors
  call Backup % Save_Int(Comm, d, vc, 'n_proc', N_Procs())

  !-----------------------------------------------------!
  !   Bulk flows and pressure drops in each direction   !
  !-----------------------------------------------------!
  call Backup % Save_Real(Comm, d, vc, 'bulk_flux_x',   bulk % flux_x)
  call Backup % Save_Real(Comm, d, vc, 'bulk_flux_y',   bulk % flux_y)
  call Backup % Save_Real(Comm, d, vc, 'bulk_flux_z',   bulk % flux_z)
  call Backup % Save_Real(Comm, d, vc, 'bulk_flux_x_o', bulk % flux_x_o)
  call Backup % Save_Real(Comm, d, vc, 'bulk_flux_y_o', bulk % flux_y_o)
  call Backup % Save_Real(Comm, d, vc, 'bulk_flux_z_o', bulk % flux_z_o)
  call Backup % Save_Real(Comm, d, vc, 'bulk_p_drop_x', bulk % p_drop_x)
  call Backup % Save_Real(Comm, d, vc, 'bulk_p_drop_y', bulk % p_drop_y)
  call Backup % Save_Real(Comm, d, vc, 'bulk_p_drop_z', bulk % p_drop_z)

  !----------------------------!
  !                            !
  !   Navier-Stokes equation   !
  !                            !
  !----------------------------!

  !--------------!
  !   Velocity   !
  !--------------!
  call Backup % Save_Variable(d, vc, 'u_velocity', Flow % u)
  call Backup % Save_Variable(d, vc, 'v_velocity', Flow % v)
  call Backup % Save_Variable(d, vc, 'w_velocity', Flow % w)

  !------------------------------------------------------!
  !   Pressure, its gradients, and pressure correction   !
  !------------------------------------------------------!
  call Backup % Save_Cell_Real(Grid, d, vc, 'press',      Flow %  p % n)
  call Backup % Save_Cell_Real(Grid, d, vc, 'press_x',    Flow %  p % x)
  call Backup % Save_Cell_Real(Grid, d, vc, 'press_y',    Flow %  p % y)
  call Backup % Save_Cell_Real(Grid, d, vc, 'press_z',    Flow %  p % z)
  call Backup % Save_Cell_Real(Grid, d, vc, 'press_corr', Flow % pp % n)

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
  call Backup % Save_Log(Comm, d, vc, 'has_pressure', Flow % has_pressure)

  !--------------!
  !              !
  !   Enthalpy   !
  !              !
  !--------------!
  if(Flow % heat_transfer) then
    call Backup % Save_Variable(d, vc, 'temp', Flow % t)
  end if

  !--------------!
  !              !
  !  Multiphase  !
  !              !
  !--------------!
  if(Flow % with_interface) then
    call Backup % Save_Variable(d, vc, 'vof_fun', Vof % fun)
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
    call Backup % Save_Variable(d, vc, 'kin', Turb % kin)
    call Backup % Save_Variable(d, vc, 'eps', Turb % eps)

    ! Other turbulent quantities
    call Backup % Save_Cell_Real(Grid, d, vc, 'p_kin',  Turb % p_kin )
    call Backup % Save_Cell_Real(Grid, d, vc, 'y_plus', Turb % y_plus)
    call Backup % Save_Cell_Real(Grid, d, vc, 'vis_t',  Turb % vis_t )
    call Backup % Save_Cell_Real(Grid, d, vc, 'vis_w',  Turb % vis_w )

    ! Turbulence quantities connected with heat transfer
    if(Flow % heat_transfer) then
      call Backup % Save_Variable(d, vc, 't2',    Turb % t2)
      call Backup % Save_Cell_Real(Grid, d, vc, 'p_t2',  Turb % p_t2 )
      call Backup % Save_Cell_Real(Grid, d, vc, 'con_w', Turb % con_w)
    end if

  end if

  !------------------------!
  !   K-eps-zeta-f model   !
  !------------------------!
  if(Turb % model .eq. K_EPS_ZETA_F .or.  &
     Turb % model .eq. HYBRID_LES_RANS) then

    ! K, eps, zeta and f22
    call Backup % Save_Variable(d, vc, 'kin',  Turb % kin)
    call Backup % Save_Variable(d, vc, 'eps',  Turb % eps)
    call Backup % Save_Variable(d, vc, 'zeta', Turb % zeta)
    call Backup % Save_Variable(d, vc, 'f22',  Turb % f22)

    ! Other turbulent quantities
    call Backup % Save_Cell_Real(Grid, d, vc,'p_kin',   Turb % p_kin  )
    call Backup % Save_Cell_Real(Grid, d, vc,'y_plus',  Turb % y_plus )
    call Backup % Save_Cell_Real(Grid, d, vc,'vis_t',   Turb % vis_t  )
    call Backup % Save_Cell_Real(Grid, d, vc,'vis_w',   Turb % vis_w  )
    call Backup % Save_Cell_Real(Grid, d, vc,'t_scale', Turb % t_scale)
    call Backup % Save_Cell_Real(Grid, d, vc,'l_scale', Turb % l_scale)

    if(Flow % heat_transfer) then
      call Backup % Save_Variable(d, vc, 't2',    Turb % t2)
      call Backup % Save_Cell_Real(Grid, d, vc, 'p_t2',  Turb % p_t2 )
      call Backup % Save_Cell_Real(Grid, d, vc, 'con_w', Turb % con_w)
    end if

  end if

  !-----------------!
  !   S-A model     !
  !-----------------!
  if(Turb % model .eq. SPALART_ALLMARAS) then

    ! K and epsilon
    call Backup % Save_Variable(d, vc, 'vis', Turb % vis)

    ! Other turbulent quantities
    call Backup % Save_Cell_Real(Grid, d, vc, 'y_plus', Turb % y_plus)
    call Backup % Save_Cell_Real(Grid, d, vc, 'vis_t',  Turb % vis_t )
    call Backup % Save_Cell_Real(Grid, d, vc, 'vis_w',  Turb % vis_w )

    ! Turbulence quantities connected with heat transfer
    if(Flow % heat_transfer) then
      call Backup % Save_Cell_Real(Grid, d, vc, 'con_w', Turb % con_w)
    end if

  end if

  !--------------!
  !   Roughness  !
  !--------------!
  call Backup % Save_Cell_Real(Grid, d, vc, 'z_o', Turb % z_o)

  !------------------!
  !   Save scalars   !
  !------------------!
  do sc = 1, Flow % n_scalars
    phi => Flow % scalar(sc)
    call Backup % Save_Variable(d, vc, phi % name, phi)
  end do

  !-----------------------------------------!
  !                                         !
  !   Turbulent statistics for all models   !
  !                                         !
  !-----------------------------------------!
  if(Turb % statistics) then

    call Backup % Save_Cell_Real(Grid, d, vc, 'u_mean', Turb % u_mean)
    call Backup % Save_Cell_Real(Grid, d, vc, 'v_mean', Turb % v_mean)
    call Backup % Save_Cell_Real(Grid, d, vc, 'w_mean', Turb % w_mean)
    call Backup % Save_Cell_Real(Grid, d, vc, 'p_mean', Turb % p_mean)
    if(Flow % heat_transfer) then
      call Backup % Save_Cell_Real(Grid, d, vc, 't_mean', Turb % t_mean)
      call Backup % Save_Cell_Real(Grid, d, vc, 'q_mean', Turb % q_mean)
    end if

    ! K and epsilon
    if(Turb % model .eq. K_EPS) then
      call Backup % Save_Cell_Real(Grid, d, vc, 'kin_mean', Turb % kin_mean)
      call Backup % Save_Cell_Real(Grid, d, vc, 'eps_mean', Turb % eps_mean)
      if(Flow % heat_transfer) then
        call Backup % Save_Cell_Real(Grid, d, vc, 't2_mean', Turb % t2_mean)
        call Backup % Save_Cell_Real(Grid, d, vc, 'ut_mean', Turb % ut_mean)
        call Backup % Save_Cell_Real(Grid, d, vc, 'vt_mean', Turb % vt_mean)
        call Backup % Save_Cell_Real(Grid, d, vc, 'wt_mean', Turb % wt_mean)
      end if
    end if

    ! K-eps-zeta-f and the hybrid model
    if(Turb % model .eq. K_EPS_ZETA_F .or.  &
       Turb % model .eq. HYBRID_LES_RANS) then
      call Backup % Save_Cell_Real(Grid, d, vc, 'kin_mean',  Turb % kin_mean )
      call Backup % Save_Cell_Real(Grid, d, vc, 'eps_mean',  Turb % eps_mean )
      call Backup % Save_Cell_Real(Grid, d, vc, 'zeta_mean', Turb % zeta_mean)
      call Backup % Save_Cell_Real(Grid, d, vc, 'f22_mean',  Turb % f22_mean )
      if(Flow % heat_transfer) then
        call Backup % Save_Cell_Real(Grid, d, vc, 't2_mean',  Turb % t2_mean)
        call Backup % Save_Cell_Real(Grid, d, vc, 'ut_mean',  Turb % ut_mean)
        call Backup % Save_Cell_Real(Grid, d, vc, 'vt_mean',  Turb % vt_mean)
        call Backup % Save_Cell_Real(Grid, d, vc, 'wt_mean',  Turb % wt_mean)
      end if
    end if

    call Backup % Save_Cell_Real(Grid, d, vc, 'uu_res', Turb % uu_res)
    call Backup % Save_Cell_Real(Grid, d, vc, 'vv_res', Turb % vv_res)
    call Backup % Save_Cell_Real(Grid, d, vc, 'ww_res', Turb % ww_res)
    call Backup % Save_Cell_Real(Grid, d, vc, 'uv_res', Turb % uv_res)
    call Backup % Save_Cell_Real(Grid, d, vc, 'uw_res', Turb % uw_res)
    call Backup % Save_Cell_Real(Grid, d, vc, 'vw_res', Turb % vw_res)

    if(Flow % heat_transfer) then
      call Backup % Save_Cell_Real(Grid, d, vc, 't2_res', Turb % t2_res)
      call Backup % Save_Cell_Real(Grid, d, vc, 'ut_res', Turb % ut_res)
      call Backup % Save_Cell_Real(Grid, d, vc, 'vt_res', Turb % vt_res)
      call Backup % Save_Cell_Real(Grid, d, vc, 'wt_res', Turb % wt_res)
    end if

    ! Scalars
    do sc = 1, Flow % n_scalars
      phi => Flow % scalar(sc)
      name_mean = phi % name
      name_mean(5:9) = '_mean'
      call Backup % Save_Cell_Real(Grid, d, vc, name_mean,  &
                                   Turb % scalar_mean(sc, :))
    end do

  end if

  !--------------------------!
  !                          !
  !   Swarm (of particles)   !
  !                          !
  !--------------------------!
  if(Flow % with_particles) then
    call Backup % Save_Swarm(Swarm, d, vc)
    call Backup % Save_Cell_Real(Grid,d,vc, 'n_deposited', Swarm % n_deposited)
    call Backup % Save_Cell_Real(Grid,d,vc, 'n_reflected', Swarm % n_reflected)
  end if

  ! Variable count (store +1 to count its own self)
  call Backup % Save_Int(Comm, d, vc, 'variable_count', vc + 1)

  if(First_Proc()) then
    print *, '# Wrote ', vc, ' variables!'
  end if

  ! Close backup file
  call Comm % Close_File(fh)

  call Profiler % Stop('Backup_Mod_Save')

  end subroutine

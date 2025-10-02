!==============================================================================!
  subroutine Save(Backup, Grid, Flow, Turb, dom)
!------------------------------------------------------------------------------!
!   Saves backup files name.backup                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Backup_Type)        :: Backup
  type(Grid_Type),   target :: Grid
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: Turb
  integer,         optional :: dom
!-----------------------------------[Locals]-----------------------------------!
  type(Comm_Type), pointer :: Comm
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: phi
  character(SL)            :: name_out, name_mean
  integer                  :: vc, sc
  integer(DP)              :: d
!==============================================================================!

  call Profiler % Start('Backup_Mod_Save')

  ! Take aliases
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
  !   (This section is very different in CPU version)   !
  !-----------------------------------------------------!
  call Backup % Save_Real(Comm, d, vc, 'bulk_u',        bulk % u)
  call Backup % Save_Real(Comm, d, vc, 'bulk_v',        bulk % v)
  call Backup % Save_Real(Comm, d, vc, 'bulk_w',        bulk % w)
  call Backup % Save_Real(Comm, d, vc, 'bulk_u_o',      bulk % u_o)
  call Backup % Save_Real(Comm, d, vc, 'bulk_v_o',      bulk % v_o)
  call Backup % Save_Real(Comm, d, vc, 'bulk_w_o',      bulk % w_o)
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
  call Backup % Save_Cell_Real(Grid, d, vc, 'press',      Flow % p  % n)
  call Backup % Save_Cell_Real(Grid, d, vc, 'press_corr', Flow % pp % n)

  !-------------------!
  !   Volume fluxes   ! -> don't use for the time being, too much trouble
  !-------------------!

  !-------------------------!
  !   Does it have outlet   !
  !-------------------------!

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

  !-----------------------!
  !                       !
  !   Turbulence models   !
  !                       !
  !-----------------------!

  !--------------!
  !   Roughness  !
  !--------------!

  !------------------!
  !   Save scalars   !
  !------------------!
  do sc = 1, Flow % n_scalars
    phi => Flow % scalar(sc)
    call Backup % Save_Variable(d, vc, phi % name, phi)
  end do

  !--------------------------!
  !                          !
  !   Swarm (of particles)   !
  !                          !
  !--------------------------!

  ! Variable count (store +1 to count its own self)
  call Backup % Save_Int(Comm, d, vc, 'variable_count', vc + 1)

  if(First_Proc()) then
    print *, '# Wrote ', vc, ' variables!'
  end if

  ! Close backup file
  call Comm % Close_File(fh)

  call Profiler % Stop('Backup_Mod_Save')

  end subroutine

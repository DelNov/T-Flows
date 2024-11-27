!==============================================================================!
  subroutine Load(Backup, Grid, Flow, Turb, bckp)
!------------------------------------------------------------------------------!
!   Loads backup files name.backup                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Backup_Type)       :: Backup
  type(Grid_Type),  target :: Grid
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  logical                  :: bckp
!-----------------------------------[Locals]-----------------------------------!
  type(Comm_Type), pointer :: Comm
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: phi
  character(SL)            :: name_in, answer, name_mean
  integer                  :: vc, sc, ts
  integer(DP)              :: d
  logical                  :: present
  real                     :: st  ! saved time, simulation time
!==============================================================================!

  ! Take aliases
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
  call Time % Set_Curr_Dt(ts)

  ! Simulation time
  call Backup % Load_Real(Comm, d, vc, 'time', st)
  call Time % Set_Time(st)

  !-----------------------------------------------------!
  !   Bulk flows and pressure drops in each direction   !
  !   (This section is very different in CPU version)   !
  !-----------------------------------------------------!
  call Backup % Load_Real(Comm, d, vc, 'bulk_u',        bulk % u)
  call Backup % Load_Real(Comm, d, vc, 'bulk_v',        bulk % v)
  call Backup % Load_Real(Comm, d, vc, 'bulk_w',        bulk % w)
  call Backup % Load_Real(Comm, d, vc, 'bulk_u_o',      bulk % u_o)
  call Backup % Load_Real(Comm, d, vc, 'bulk_v_o',      bulk % v_o)
  call Backup % Load_Real(Comm, d, vc, 'bulk_w_o',      bulk % w_o)
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

  !--------------------------------------!
  !   Pressure and pressure correction   !
  !--------------------------------------!
  call Backup % Load_Cell_Real(Grid, d, vc, 'press',      Flow % p  % n)
  call Backup % Load_Cell_Real(Grid, d, vc, 'press_corr', Flow % pp % n)

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
    call Backup % Load_Variable(d, vc, 'temp', Flow, Flow % t)
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
  !   Load scalars   !
  !------------------!
  do sc = 1, Flow % n_scalars
    phi => Flow % scalar(sc)
    call Backup % Load_Variable(d, vc, phi % name, Flow, phi)
  end do

  !--------------------------!
  !                          !
  !   Swarm (of particles)   !
  !                          !
  !--------------------------!

  ! Close backup file
  call Comm % Close_File(fh)

  call Profiler % Stop('Backup_Mod_Load')

  end subroutine

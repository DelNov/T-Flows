!==============================================================================!
  program Processor
!------------------------------------------------------------------------------!
!   Unstructured finite volume 'LES'/RANS solver.                              !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Const_Mod
  use Comm_Mod
  use Cpu_Timer_Mod
  use File_Mod,      only: problem_name
  use Field_Mod,     only: Field_Type, Field_Mod_Allocate, heat_transfer
  use Turb_Mod
  use Grid_Mod
  use Eddies_Mod
  use Bulk_Mod
  use Var_Mod,       only: Var_Type
  use Solver_Mod,    only: Solver_Mod_Create
  use Info_Mod
  use Work_Mod,      only: Work_Mod_Allocate
  use User_Mod
  use Save_Results_Mod
  use Control_Mod
  use Monitor_Mod
  use Backup_Mod
  use Surf_Mod
  use Multiphase_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Calling]------------------------------------!
  real :: Correct_Velocity
!----------------------------------[Locals]------------------------------------!
  integer               :: n, sc, tp
  real                  :: mass_res
  logical               :: backup, save_now, exit_now
  type(Grid_Type)       :: grid            ! grid used in computations
  type(Field_Type)      :: flow            ! flow field we will be solving for
  type(Swarm_Type)      :: swarm           ! swarm of particles
  type(Surf_Type)       :: surf            ! interface between two phases
  type(Turb_Type)       :: turb            ! turbulence modelling
  type(Multiphase_Type) :: mult            ! multiphase modelling
  type(Solver_Type)     :: sol             ! linear solvers
  type(Turb_Plane_Type) :: turb_planes     ! holder for synthetic turbulences
  real                  :: time            ! physical time of the simulation
  integer               :: first_dt        ! first time step in this run
  integer               :: last_dt         ! number of time steps
  integer               :: max_ini         ! max number of inner iterations
  integer               :: min_ini         ! min number of inner iterations
  integer               :: n_stat          ! starting time step for statistic
  integer               :: ini             ! inner iteration counter
  integer               :: bsi, rsi        ! backup and results save interval
  real                  :: simple_tol      ! tolerance for SIMPLE algorithm
  integer               :: sc_cr           ! system clock count rate
  integer               :: sc_ini, sc_cur  ! system clock start and end rate
  real                  :: wt_max
!==============================================================================!

  ! Initialize program profler
  call Cpu_Timer_Mod_Start('Main')

  ! Get starting time
  call system_clock(count_rate=sc_cr)  ! get system clock clock rate
  call system_clock(sc_ini)            ! get system clock initial count rate

  ! Initialize variables
  time =  0.
  backup = .false.  ! can turn .true. in Load_Backup

  !------------------------------!
  !   Start parallel execution   !
  !------------------------------!
  call Comm_Mod_Start

  !--------------------------------!
  !   Splash out the logo screen   !
  !--------------------------------!
  if(this_proc  < 2) then
    call Logo_Pro
  end if

  !---------------------------------------------!
  !   Open control file and read problem name   !
  !---------------------------------------------!
  call Control_Mod_Open_File('control')
  call Control_Mod_Problem_Name(problem_name)

  ! Load the finite volume grid
  call Grid_Mod_Load_Cns(grid, this_proc)

  ! Allocate memory for working arrays and comm.
  call Work_Mod_Allocate(grid, rc=30, rf=1, rn=1, ic=4, if=0, in=1)
  call Grid_Mod_Allocate_Comm(grid)

  call Grid_Mod_Load_Geo(grid, this_proc)
  call Grid_Mod_Create_Buffers(grid)
  call Grid_Mod_Load_Maps(grid)       ! maps should move to .cns file soon

  call Comm_Mod_Wait

  ! Get the number of time steps from the control file
  call Control_Mod_Number_Of_Time_Steps(last_dt, verbose=.true.)
  call Control_Mod_Starting_Time_Step_For_Statistics(n_stat, verbose=.true.)

  ! Read physical models from control file
  call Read_Control_Physical(flow, turb, mult, swarm)

  ! Allocate memory for all variables
  call Field_Mod_Allocate(flow, grid)
  call Turb_Mod_Allocate(turb, flow)
  call Swarm_Mod_Allocate(swarm, flow)
  call Multiphase_Mod_Allocate(mult, flow)
  call User_Mod_Allocate(grid)

  ! Read numerical models from control file (after the memory is allocated)
  call Read_Control_Numerical(flow, turb, mult)

  call Grid_Mod_Calculate_Face_Geometry(grid)
  call Grid_Mod_Find_Nodes_Cells(grid)         ! for Lagrangian particle track
  call Grid_Mod_Find_Periodic_Faces(grid)      ! for Lagrangian particle track
  call Grid_Mod_Find_Cells_Faces(grid)         ! for Multiphase Module
  call Grid_Mod_Calculate_Global_Volumes(grid)

  ! Allocate memory for linear systems of equations
  ! (You need face geomtry for this step)
  call Solver_Mod_Create(sol, grid)

  call Load_Physical_Properties(flow, mult, swarm)

  call Load_Boundary_Conditions(flow, turb, mult, turb_planes)

  ! First time step is one, unless read from backup otherwise
  first_dt = 0

  ! Read backup file if directed so, and set the "backup" to .true. or .false.
  call Backup_Mod_Load(flow, swarm, turb, mult, first_dt, n_stat, backup) 

  ! Initialize variables
  if(.not. backup) then
    call Initialize_Variables(flow, turb, mult, swarm)
    if(multiphase_model .eq. VOLUME_OF_FLUID) then
      call Multiphase_Mod_Update_Physical_Properties(mult)
    end if
  end if

  ! Initialize monitoring points
  call Monitor_Mod_Initialize(grid, backup)

  ! Plane for calcution of overall mass fluxes
  call Control_Mod_Point_For_Monitoring_Planes(flow % bulk % xp,  &
                                               flow % bulk % yp,  &
                                               flow % bulk % zp)

  ! Prepare ...
  call Bulk_Mod_Monitoring_Planes_Areas(flow % bulk, grid)

  if(turbulence_model .eq. LES_SMAGORINSKY .and. .not. backup) then
    call Find_Nearest_Wall_Cell(turb)
  end if

  if(turbulence_model .eq. HYBRID_LES_PRANDTL .and. .not. backup) then
    call Find_Nearest_Wall_Cell(turb)
  end if

  ! Prepare the gradient matrix for velocities
  call Field_Mod_Calculate_Grad_Matrix(flow)

  ! Print the areas of monitoring planes
  call Bulk_Mod_Print_Areas(flow % bulk)

  ! Compute deltas for Spalart-Allmaras models
  call Turb_Mod_Calculate_Deltas(turb)

  !---------------!
  !               !
  !   Time loop   !
  !               !
  !---------------!

  call Control_Mod_Time_Step(flow % dt, verbose=.true.)
  call Control_Mod_Backup_Save_Interval (bsi,    verbose=.true.)
  call Control_Mod_Results_Save_Interval(rsi,    verbose=.true.)
  call Control_Mod_Wall_Time_Max_Hours  (wt_max, verbose=.true.)
  wt_max = wt_max * 3600  ! make it in seconds

  ! It will save results in .vtk or .cgns file format,
  ! depending on how the code was compiled
  call Save_Results(flow, turb, mult, swarm, first_dt, .true.)   ! save inside
  call Save_Results(flow, turb, mult, swarm, first_dt, .false.)  ! save boundary
  call Save_Swarm(swarm, first_dt)

  do n = first_dt + 1, last_dt

    ! Update turbulent planes
    do tp = 1, turb_planes % n_planes
      call Eddies_Mod_Superimpose(turb_planes % plane(tp))
      call Eddies_Mod_Advance    (turb_planes % plane(tp))
    end do

    time = time + flow % dt

    ! Beginning of time step
    call User_Mod_Beginning_Of_Time_Step(flow, turb, mult, swarm, n, time)

    ! Start info boxes.
    call Info_Mod_Time_Start()
    call Info_Mod_Iter_Start()
    call Info_Mod_Bulk_Start()

    ! Initialize and print time info box
    call system_clock(sc_cur)
    call Info_Mod_Time_Fill( n, time, real(sc_cur-sc_ini)/real(sc_cr) )
    call Info_Mod_Time_Print()

    ! Turbulence models initializations
    call Turb_Mod_Init(turb, sol)

    !--------------------------!
    !   Inner-iteration loop   !
    !--------------------------!
    call Control_Mod_Max_Simple_Iterations(max_ini)
    call Control_Mod_Min_Simple_Iterations(min_ini)

    ! Volume of Fluid
    if(multiphase_model .eq. VOLUME_OF_FLUID) then
      flow % m_flux % o = flow % m_flux % n / flow % density_f
      ! Update the values at boundaries
      call Update_Boundary_Values(flow, turb, mult)
      call Multiphase_Mod_Compute_Vof(mult, sol, flow % dt, n)
    else
      flow % m_flux % o = flow % m_flux % n
    end if

    do ini = 1, max_ini

      ! Beginning of iteration
      call User_Mod_Beginning_Of_Iteration(flow, turb, mult, swarm, n, time)

      call Info_Mod_Iter_Fill(ini)

      call Field_Mod_Grad_Pressure(flow, flow % p,  &
                                   flow % density,  &
                                   grav_x, grav_y, grav_z)

      ! Compute velocity gradients
      call Field_Mod_Grad_Variable(flow, flow % u)
      call Field_Mod_Grad_Variable(flow, flow % v)
      call Field_Mod_Grad_Variable(flow, flow % w)

      ! All velocity components one after another
      call Compute_Momentum(flow, turb, mult, 1, sol, flow % dt, ini)
      call Compute_Momentum(flow, turb, mult, 2, sol, flow % dt, ini)
      call Compute_Momentum(flow, turb, mult, 3, sol, flow % dt, ini)

      ! Refresh buffers for a % sav before discretizing for pressure
      ! (Can this call be somewhere in Compute Pressure?)
      call Grid_Mod_Exchange_Real(grid, sol % a % sav)

      call Balance_Mass(flow)
      call Compute_Pressure(flow, mult, sol, flow % dt, ini)

      call Field_Mod_Grad_Pressure_Correction(flow, flow % pp)

      call Field_Mod_Calculate_Fluxes(flow, flow % m_flux % n)
      mass_res = Correct_Velocity(flow, sol, flow % dt, ini)

      ! Energy (practically temperature)
      if(heat_transfer) then
        call Compute_Energy(flow, turb, mult, sol, flow % dt, ini)
      end if

      ! Passive scalars
      do sc = 1, flow % n_scalars
        call Compute_Scalar(flow, turb, mult, sol, flow % dt, ini, sc)
      end do

      ! Deal with turbulence (if you dare ;-))
      call Turb_Mod_Main(turb, sol, n, ini)

      ! Update the values at boundaries
      call Convective_Outflow(flow, turb, mult, flow % dt)
      call Update_Boundary_Values(flow, turb, mult)

      ! End of the current iteration
      call Info_Mod_Iter_Print()

      ! End of iteration
      call User_Mod_End_Of_Iteration(flow, turb, mult, swarm, n, time)

      if(ini >= min_ini) then
        call Control_Mod_Tolerance_For_Simple_Algorithm(simple_tol)
        if( flow % u  % res <= simple_tol .and.  &
            flow % v  % res <= simple_tol .and.  &
            flow % w  % res <= simple_tol .and.  &
            mass_res        <= simple_tol ) goto 1
      end if

    end do

    ! End of the current time step
1   call Info_Mod_Bulk_Print()

    ! Write the values in monitoring points
    if(.not. heat_transfer) then
      call Monitor_Mod_Write_4_Vars(n, flow)
    else
      call Monitor_Mod_Write_5_Vars(n, flow)
    end if

    ! Calculate mean values
    call Turb_Mod_Calculate_Mean(turb, n_stat, n)
    call User_Mod_Calculate_Mean(turb, n_stat, n)

    ! Adjust pressure drops to keep the mass fluxes constant
    call Bulk_Mod_Adjust_P_Drops(flow % bulk, flow % dt)

    ! Just before the end of time step
    call User_Mod_End_Of_Time_Step(flow, turb, mult, swarm, n, time)

    !----------------------!
    !   Save the results   !
    !----------------------!
    inquire(file='exit_now', exist=exit_now)
    inquire(file='save_now', exist=save_now)

    ! Is it time to save the backup file?
    if(save_now           .or.  &
       exit_now           .or.  &
       mod(n, bsi) .eq. 0 .or.  &
       real(sc_cur-sc_ini)/real(sc_cr) > wt_max) then
      call Backup_Mod_Save(flow, swarm, turb, mult, n, n_stat)
    end if

    ! Is it time to save results for post-processing
    if(save_now           .or.  &
       exit_now           .or.  &
       mod(n, rsi) .eq. 0 .or.  &
       real(sc_cur-sc_ini)/real(sc_cr) > wt_max) then
      call Comm_Mod_Wait
      call Save_Results(flow, turb, mult, swarm, n, .true.)   ! save inside
      call Save_Results(flow, turb, mult, swarm, n, .false.)  ! save bnd
      call Save_Swarm(swarm, n)

      if(multiphase_model .eq. VOLUME_OF_FLUID) then
!        call Surf_Mod_Allocate(surf, flow)
!        call Surf_Mod_Place_At_Var_Value(surf,        &
!                                         mult % vof,  &
!                                         sol,         &
!                                         0.5,         &
!                                         .false.)  ! don't print messages
!        call Save_Surf(surf, n)
!        call Surf_Mod_Clean(surf)
      end if

      ! Write results in user-customized format
      call User_Mod_Save_Results(flow, turb, mult, swarm, n)
      ! call User_Mod_Save_Swarm(swarm, n)  to be done!
    end if

    if(save_now) then
      if(this_proc < 2) then
        open (9, file='save_now', status='old')
        close(9, status='delete')
      end if
    end if

    if(exit_now) then
      if(this_proc < 2) then
        open (9, file='exit_now', status='old')
        close(9, status='delete')
      end if
      goto 2
    end if

    ! Ran more than a set limit
    if(real(sc_cur-sc_ini)/real(sc_cr) > wt_max) then
      goto 2
    end if

  end do ! n, number of time steps

  ! After the time loop, decrease "n" since ...
  ! ... it is one above the loop boundaries here
  n = n - 1

  ! User function for end of simulation
  call User_Mod_End_Of_Simulation(flow, turb, mult, swarm, n, time)

  ! Save backup and post-processing files at exit
  call Comm_Mod_Wait
  call Backup_Mod_Save(flow, swarm, turb, mult, n, n_stat)
  call Save_Results(flow, turb, mult, swarm, n, .true.)   ! save inside
  call Save_Results(flow, turb, mult, swarm, n, .false.)  ! save bnd

  ! Write results in user-customized format
  call User_Mod_Save_Results(flow, turb, mult, swarm, n)

  if(this_proc < 2) then
    open(9, file='stop')
    close(9)
  end if

2 if(this_proc  < 2) print *, '# Exiting !'

  ! Close monitoring files
  call Monitor_Mod_Finalize

  ! Make the final call to user function
  call User_Mod_Before_Exit(grid)

  call Cpu_Timer_Mod_Stop('Main')
  call Cpu_Timer_Mod_Statistics

  !----------------------------!
  !   End parallel execution   !
  !----------------------------!
  call Comm_Mod_End

  end program

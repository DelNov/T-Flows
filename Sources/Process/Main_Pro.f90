#include "../Shared/Assert.h90"

!==============================================================================!
  program Process_Prog
!------------------------------------------------------------------------------!
!>  This is the main function of the T-Flows' sub-program Process.  It manages
!>  the entire workflow of CFD simulations. It processes the complex interplay
!>  of initializing the computing environment, preparing grids, initializing
!>  solvers, handling domain interactions, and managing data across various
!>  computational stages. From executing core computational routines to
!>  conducting post-processing and smoothly transitioning to program exit,
!>  this function stands as a critical component of Process, ensuring that
!>  every step of the simulation is executed correctly and in synchronization
!>  with other parts of the code.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization:                                                          !
!     - Sets up profiler for performance tracking.                             !
!     - Initializes control file names for the root and individual domains.    !
!     - Initializes time and backup-related variables.                         !
!     - Starts parallel execution if necessary.                                !
!     - Displays the program logo.                                             !
!   * Control file processing:                                                 !
!     - Opens and reads control files for root and domain-specific settings.   !
!     - Initializes Info_Mod for managing information and logging.             !
!     Grid and field preparation:                                              !
!     - Loads and prepares grids for each domain.                              !
!     - Allocates memory for working arrays, essential for RSM models.         !
!     - Sets up time steps and reads control information for physical models.  !
!   * Memory allocation for variables:                                         !
!     - Allocates memory for flow fields, turbulence, VOF, swarms, and solvers.!
!     - Initializes numerical schemes and boundary conditions.                 !
!   * Interface creation and backup loading:                                   !
!     - Creates interfaces between domains.                                    !
!     - Loads backup files if directed and initializes variables accordingly.  !
!   * Time loop setup:                                                         !
!     - Sets intervals for backup and results saving.                          !
!     - Computes wall distance and performs potential initialization.          !
!   * Main time loop:                                                          !
!     - Iteratively processes each time step.                                  !
!     - Handles turbulence models, interface tracking, and                     !
!       Lagrangian particle tracking model.                                    !
!     - Executes the PISO algorithm for pressure-velocity coupling.            !
!     - Updates boundary values and solves for energy and scalar fields.       !
!     - Manages data exchange between domains.                                 !
!     - Performs post-iteration tasks, including user-defined functions.       !
!   * Result saving and exiting:                                               !
!     - Saves results at specified intervals.                                  !
!     - Checks for exit conditions based on time or external signals.          !
!     - Calls user functions at the end of the simulation.                     !
!   * Finalization:                                                            !
!     - Closes monitoring files and performs final user function calls.        !
!     - Finalizes work arrays, profiler and solvers.                           !
!     - Ends parallel execution and exits the program.                         !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Process_Mod
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Locals]------------------------------------!
  character(len=7)      :: root_control    = 'control'
  character(len=9)      :: dom_control(MD) = 'control.d'
  integer               :: sc, tp, ldt
  logical               :: read_backup(MD), exit_now, pot_init
  type(Grid_Type)       :: Grid(MD)        ! grid used in computations
  type(Field_Type)      :: Flow(MD)        ! flow field we will be solving for
  type(Swarm_Type)      :: Swarm(MD)       ! swarm of particles
  type(Turb_Type)       :: Turb(MD)        ! turbulence modelling
  type(Vof_Type)        :: Vof(MD)         ! multiphase modelling with vof
  type(Solver_Type)     :: Sol(MD)         ! native and PETSc linear solvers
  type(Turb_Plane_Type) :: turb_planes(MD) ! holder for synthetic turbulences
  type(Monitor_Type)    :: monitor(MD)     ! monitors
  type(Porosity_Type)   :: Por(MD)         ! porosity
  type(Interface_Type)  :: inter(MD,MD)    ! interfaces between domains
  integer               :: n_stat_t        ! first time step for turb. statistic
  integer               :: n_stat_p        ! first time step for swarm statistic
  integer               :: first_dt_p      ! first t.s. for swarm computation
  integer               :: n_dom, d        ! number of domains, domain counter
!==============================================================================!

  ! Initialize program profler
  call Profiler % Start('Main')

  ! Initialize control file names
  root_control = 'control'             ! root control file name
  do d = 1, MD                         ! domain control file names
    write(dom_control(d), '(a8,i1)') 'control.', d
  end do

  ! Initialize variables
  call Time % Set_Time(0.0)  ! initialize time to zero
  read_backup(:) = .false.   ! can turn .true. in Load_Backup

  !------------------------------!
  !   Start parallel execution   !
  !------------------------------!
  call Global % Start_Parallel

  !--------------------------------!
  !   Splash out the logo screen   !
  !--------------------------------!
  call Process % Logo_Pro()

  !-----------------------!
  !   Open control file   !
  !-----------------------!
  call Control % Open_Root_File(root_control)

  call Control % Number_Of_Domains(n_dom)
  if(n_dom > 1) then
    do d = 1, n_dom
      call Control % Open_Domain_File(d, dom_control(d))
    end do
  end if

  !-------------------------!
  !   Initialize Info_Mod   !
  !-------------------------!
  call Info % Start_Info()

  !-----------------------------------------------!
  !   Load and prepare all grids for processing   !
  !-----------------------------------------------!
  do d = 1, n_dom
    call Control % Switch_To_Domain(d)  ! take domain's d control file
    call Grid(d) % Load_And_Prepare_For_Processing(d)
  end do
  call Control % Switch_To_Root()  ! out of domain loop - go back to root

  ! Allocate memory for working arrays (RSM models are memory hungry)
  call Work % Allocate_Work(Grid, n_r_cell=24,  n_r_face=8,  n_r_node=8,  &
                                  n_i_cell= 8,  n_i_face=8,  n_i_node=8)

  ! Initialize first and current and read the last time step
  call Time % Set_Curr_Dt(0)
  call Time % Set_First_Dt(0)
  call Control % Number_Of_Time_Steps(ldt, verbose=.true.)
  call Time % Set_Last_Dt(ldt)
  call Control % Starting_Time_Step_For_Turb_Statistics(n_stat_t,  &
                                                        verbose = .true.)
  call Control % Starting_Time_Step_For_Swarm_Statistics(n_stat_p,  &
                                                         verbose = .true.)
  call Control % Starting_Time_Step_For_Swarm_Computation(first_dt_p,  &
                                                          verbose = .true.)

  ! Read physical models for each domain from control file
  do d = 1, n_dom
    call Control % Switch_To_Domain(d)  ! take proper control file
    call Read_Control % Physical_Models(Flow(d), Turb(d), Vof(d), Swarm(d))
  end do

  !----------------------------------------------------------!
  !   Allocate memory for all variables (over all domains)   !
  !----------------------------------------------------------!
  do d = 1, n_dom
    call Control % Switch_To_Domain(d)  ! take proper control file
    call Read_Control % Solvers(Flow(d), Turb(d), Vof(d), Sol(d))
    call Sol(d)   % Create_Solver(Grid(d))
    call Flow(d)  % Create_Field(Sol(d) % Nat % A)  ! will store pnt_matrix
    call Turb(d)  % Create_Turb(Flow(d))
    call Vof(d)   % Create_Vof(Flow(d))
    call Swarm(d) % Create_Swarm(Flow(d), Turb(d), Vof(d))

    ! Read time step from root
    call Control % Switch_To_Root()
    call Control % Time_Step(Flow(d) % dt, verbose=.true.)
    call Control % Switch_To_Domain(d)  ! go back to local domain's control

    ! Read numerical models from control file (after the memory is allocated)
    call Read_Control % Numerical_Schemes(Flow(d), Turb(d), Vof(d))

    call Grid(d) % Find_Nodes_Cells()
    call Grid(d) % Calculate_Weights_Cells_To_Nodes()  ! needed for front
    call Grid(d) % Calculate_Global_Volumes()
    call Flow(d) % Calculate_Grad_Matrix()

    call Read_Control % Physical_Properties(Flow(d), Vof(d), Swarm(d))
    call Read_Control % Boundary_Conditions(Flow(d), Turb(d), Vof(d),   &
                                            turb_planes(d))
  end do

  ! Create interfaces
  call Control % Switch_To_Root()
  call Interface_Mod_Create(inter, Grid, n_dom)

  ! Read backup file if directed so, and set the "backup" to .true. or .false.
  do d = 1, n_dom
    call Control % Switch_To_Domain(d)  ! take proper control file
    call Backup % Load(Flow(d), Turb(d), Vof(d), Swarm(d), read_backup(d))

    ! Initialize variables
    if(.not. read_backup(d)) then
      call Process % Initialize_Variables(Flow(d), Turb(d),  &
                                          Vof(d), Swarm(d), Sol(d))
    end if

    if(Flow(d) % with_interface) then
      if (read_backup(d))  then
        Flow % inside_piso_loop = .false.
      end if
      call Vof(d) % Update_Physical_Properties()
    end if

    ! Initialize monitoring points
    call Monitor(d) % Initialize(Flow(d), read_backup(d), domain=d)

    call Por(d) % Create_Porosity(Grid(d))

    ! Plane for calcution of overall mass fluxes
    call Control % Point_For_Monitoring_Planes(Flow(d) % bulk % xp,  &
                                               Flow(d) % bulk % yp,  &
                                               Flow(d) % bulk % zp)
    ! Prepare ...
    call Bulk_Mod_Monitoring_Planes_Areas(Flow(d) % bulk, Grid(d))

    ! Print the areas of monitoring planes
    call Bulk_Mod_Print_Areas(Flow(d) % bulk)

    ! Compute deltas for Spalart-Allmaras models
    call Turb(d) % Calculate_Deltas()

  end do

  !---------------!
  !               !
  !   Time loop   !
  !               !
  !---------------!
  call Control % Switch_To_Root()
  call Control % Backup_Save_Interval  (backup % interval, verbose=.true.)
  call Control % Results_Save_Interval (Results % interval, verbose=.true.)
  call Control % Save_Initial_Condition(Results % initial,  verbose=.true.)
  call Control % Save_Results_At_Boundaries(Results % boundary)
  call Control % Swarm_Save_Interval(Results % interval_swarm, verbose=.true.)

  !-------------------------------------------------------!
  !   Compute wall distance - it is not saved in backup   !
  !      file and could be set to -1.0 in .dim file       !
  !-------------------------------------------------------!
  do d = 1, n_dom
    call Flow(d) % Compute_Wall_Distance(Sol(d))
  end do

  !-------------------------------------------------------------!
  !   Perform potential initialization in the first time step   !
  !-------------------------------------------------------------!
  if(Time % First_Dt() .eq. 0) then
    do d = 1, n_dom
      call Control % Switch_To_Domain(d)  ! not sure if this call is needed
      call Control % Potential_Initialization(pot_init, .true.)
      call Flow(d) % Potential_Initialisation(Sol(d), pot_init)
    end do
  end if

  ! Good time to call user function for beginning of simulation
  do d = 1, n_dom
    call User_Mod_Beginning_Of_Simulation(Flow(d), Turb(d), Vof(d), Swarm(d))
  end do

  ! Save initial condition
  call Results % Main_Results(n_dom, Flow, Turb, Vof, Swarm, exit_now)

  !-------------------------------------!
  !   The time loop really begins now   !
  !-------------------------------------!
  do while (Time % Needs_More_Steps())

    !------------------------------------!
    !   Preparations for new time step   !
    !------------------------------------!
    do d = 1, n_dom

      call Control % Switch_To_Domain(d)  ! not sure if this call is needed

      ! Update turbulent planes
      do tp = 1, turb_planes(d) % n_planes
        call Eddies_Mod_Superimpose(turb_planes(d) % plane(tp))
        call Eddies_Mod_Advance    (turb_planes(d) % plane(tp))
      end do

      if(d .eq. 1) call Time % Increase_Time(Flow(d) % dt)

      ! Beginning of time step
      call User_Mod_Beginning_Of_Time_Step(Flow(d), Turb(d), Vof(d), Swarm(d))

      ! Start info boxes
      call Info % Time_Start()
      call Info % Iter_Start()
      call Info % Bulk_Start()

      ! Initialize and print time info box
      if(d .eq. 1) then
        call Info % Time_Fill(Time % Curr_Dt(), Time % Get_Time())
        call Info % Time_Print()
      end if

      ! Turbulence models initializations
      call Turb(d) % Init_Turb()

      ! Interface tracking
      if(Flow(d) % with_interface) then
        call Process % Update_Boundary_Values(Flow(d), Turb(d), Vof(d), 'VOF')
        call Vof(d) % Main_Vof(Flow(d), Turb(d), Sol(d))
        call Vof(d) % Update_Physical_Properties()
      end if

      ! Lagrangian particle tracking
      if(Flow(d) % with_particles) then
        call User_Mod_Insert_Particles(Flow(d), Turb(d), Vof(d), Swarm(d))
      end if

    end do  ! through domains

    !--------------------------!
    !   Inner-iteration loop   !
    !--------------------------!
    call Control % Switch_To_Root()
    call Read_Control % Iterations()

    do while (Iter % Needs_More_Iterations(Flow, n_dom))

      ! Exchange data between domains
      call User_Mod_Interface_Exchange(inter, Flow, Turb, Vof, Swarm, n_dom)

      do d = 1, n_dom

        call Control % Switch_To_Domain(d)

        ! Beginning of iteration
        call User_Mod_Beginning_Of_Iteration(Flow(d), Turb(d), Vof(d), Swarm(d))

        call Info % Iter_Fill(Iter % Current())

        ! Future? call Process % Simple_Step(Flow(d), Turb(d), Vof(d), Sol(d))

        ! Compute velocity gradients
        call Flow(d) % Grad_Variable(Flow(d) % u)
        call Flow(d) % Grad_Variable(Flow(d) % v)
        call Flow(d) % Grad_Variable(Flow(d) % w)

        ! All three velocity components one after another
        call Process % Compute_Momentum(Flow(d), Turb(d), Vof(d), Por(d),  &
                                        Sol(d))
        call Process % Compute_Pressure(Flow(d), Vof(d), Sol(d))
        call Process % Correct_Velocity(Flow(d), Vof(d), Sol(d))

        call Process % Piso_Algorithm(Flow(d), Turb(d), Vof(d), Por(d), Sol(d))

        call Flow(d) % Calculate_Bulk_Fluxes(Flow(d) % v_flux % n)

        ! Deal with turbulence
        call Turb(d) % Main_Turb(Sol(d))

        ! Energy (practically temperature)
        call Process % Compute_Energy(Flow(d), Turb(d), Vof(d), Sol(d))

        ! Passive scalars
        do sc = 1, Flow(d) % n_scalars
          call Process % Compute_Scalar(Flow(d), Turb(d), Vof(d), Sol(d), sc)
        end do

        ! Update the values at boundaries
        call Process % Update_Boundary_Values(Flow(d), Turb(d), Vof(d), 'ALL')

        ! End of the current iteration
        call Info % Iter_Print(d)

        ! End of iteration
        call User_Mod_End_Of_Iteration(Flow(d), Turb(d), Vof(d), Swarm(d))

      end do  ! through domains

    end do    ! through outer iterations

    !----------------------------------!
    !   End of the current time step   !
    !----------------------------------!
    do d = 1, n_dom
      call Process % Convective_Outflow(Flow(d), Turb(d), Vof(d))
    end do

    do d = 1, n_dom
      call Info % Bulk_Print(Flow(d), d, n_dom)
    end do

    do d = 1, n_dom
      call Control % Switch_To_Domain(d)

      ! Write the values in monitoring points
      call Monitor(d) % Write_Vars(Flow(d), Time % Curr_Dt())

      ! Calculate mean values
      call Turb(d) % Calculate_Mean(n_stat_t)
      call User_Mod_Calculate_Mean(Turb(d), n_stat_t)

      ! Adjust pressure drops to keep the mass fluxes constant
      call Bulk_Mod_Adjust_P_Drops(Flow(d) % bulk, Flow(d) % dt)

      ! Lagrangian particle tracking
      if(Flow(d) % with_particles) then
        if(Time % Curr_Dt() >= first_dt_p) then
          call Swarm(d) % Advance_Particles(n_stat_p, first_dt_p)
        end if
      end if

      ! Just before the end of time step
      call User_Mod_End_Of_Time_Step(Flow(d), Turb(d), Vof(d), Swarm(d),  &
                                     n_stat_t, n_stat_p)
    end do

    !----------------------!
    !   Save the results   !
    !----------------------!
    call Results % Main_Results(n_dom, Flow, Turb, Vof, Swarm, exit_now)

    ! Ran more than a set wall clock time limit
    if(Info % Time_To_Exit() .or. exit_now) then
      goto 2
    end if

    ! Last time step reached; call user function for end of simulation
    if(Time % Curr_Dt() .eq. Time % Last_Dt()) then
      do d = 1, n_dom
        call Control % Switch_To_Domain(d)
        call User_Mod_End_Of_Simulation(Flow(d), Turb(d), Vof(d), Swarm(d))
      end do
    end if

  end do ! curr_dt until the last time step

2 if(First_Proc()) print *, '# Exiting !'

  do d = 1, n_dom
    ! Close monitoring files
    call Monitor(d) % Finalize()

    ! Make the final call to user function
    call User_Mod_Before_Exit(Grid(d))
  end do

  ! Looks a bit ugly
  call Work % Finalize_Work()

  ! Finalize the program profler
  call Profiler % Stop('Main')
  call Profiler % Statistics(indent=33)

  !--------------------------!
  !   Finalize the solvers   !
  !--------------------------!
  do d = 1, n_dom
    call Sol(d) % End()
  end do

  !----------------------------!
  !   End parallel execution   !
  !----------------------------!
  call Global % End_Parallel

  end program

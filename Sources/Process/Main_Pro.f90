!==============================================================================!
  program Processor
!------------------------------------------------------------------------------!
!   Unstructured finite volume 'LES'/RANS solver.                              !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Eddies_Mod
  use Work_Mod,      only: Work_Mod_Allocate
  use User_Mod
  use Results_Mod
  use Backup_Mod
  use Monitor_Mod
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Locals]------------------------------------!
  character(len=7)      :: root_control    = 'control'
  character(len=9)      :: dom_control(MD) = 'control.d'
  integer               :: curr_dt, sc, tp
  logical               :: read_backup(MD), exit_now, pot_init
  type(Grid_Type)       :: grid(MD)        ! grid used in computations
  type(Field_Type)      :: flow(MD)        ! flow field we will be solving for
  type(Swarm_Type)      :: swarm(MD)       ! swarm of particles
  type(Turb_Type)       :: turb(MD)        ! turbulence modelling
  type(Multiphase_Type) :: mult(MD)        ! multiphase modelling
  type(Solver_Type)     :: Sol(MD)         ! linear solvers
  type(Turb_Plane_Type) :: turb_planes(MD) ! holder for synthetic turbulences
  type(Monitor_Type)    :: monitor(MD)     ! monitors
  type(Interface_Type)  :: inter(MD,MD)    ! interfaces between domains
  real                  :: time            ! physical time of the simulation
  integer               :: first_dt        ! first time step in this run
  integer               :: last_dt         ! number of time steps
  integer               :: max_ini         ! max number of inner iterations
  integer               :: min_ini         ! min number of inner iterations
  integer               :: n_stat_t        ! first time step for turb. statistic
  integer               :: n_stat_p        ! first time step for swarm statistic
  integer               :: first_dt_p      ! first t.s. for swarm computation
  integer               :: ini             ! inner iteration counter
  integer               :: prsi            ! particles save interval
  real                  :: simple_tol      ! tolerance for SIMPLE algorithm
  integer               :: n_dom           ! number of domains
  integer               :: d               ! domain counter
!==============================================================================!

  ! Initialize program profler
  call Cpu_Timer % Start('Main')

  ! Initialize control file names
  root_control = 'control'             ! root control file name
  do d = 1, MD                         ! domain control file names
    write(dom_control(d), '(a8,i1)') 'control.', d
  end do

  ! Initialize variables
  time           =  0.      ! initialize time to zero
  read_backup(:) = .false.  ! can turn .true. in Load_Backup

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

  !-----------------------!
  !   Open control file   !
  !-----------------------!
  call Control_Mod_Open_Root_File(root_control)

  call Control_Mod_Number_Of_Domains(n_dom)
  if(n_dom > 1) then
    do d = 1, n_dom
      call Control_Mod_Open_Domain_File(d, dom_control(d))
    end do
  end if

  !-------------------------!
  !   Initialize Info_Mod   !
  !-------------------------!
  call Info_Mod_Start()

  !--------------------!
  !   Read all grids   !
  !--------------------!
  do d = 1, n_dom
    call Control_Mod_Switch_To_Domain(d)  ! take domain's d control file
    call Control_Mod_Problem_Name(problem_name(d))

    ! Load the finite volume grid
    call Grid_Mod_Load_Cfn(grid(d), this_proc, domain=d)
    call Grid_Mod_Load_Dim(grid(d), this_proc, domain=d)
    call Grid_Mod_Calculate_Face_Geometry(grid(d))

    call Grid_Mod_Form_Cells_Comm(grid(d))
    call Grid_Mod_Form_Maps(grid(d))

    call Comm_Mod_Wait

    call Save_Vtu_Faces(grid(d))
    call Save_Vtu_Faces(grid(d), plot_shadows=.true.)
  end do

  ! Out of domain loop - go back to root
  call Control_Mod_Switch_To_Root()

  ! Allocate memory for working arrays
  call Work_Mod_Allocate(grid, rc=30, rf=6, rn=12, ic=4, if=6, in=1)

  ! Get the number of time steps from the control file
  call Control_Mod_Number_Of_Time_Steps(last_dt, verbose=.true.)
  call Control_Mod_Starting_Time_Step_For_Turb_Statistics(n_stat_t,  &
                                                          verbose = .true.)
  call Control_Mod_Starting_Time_Step_For_Swarm_Statistics(n_stat_p,  &
                                                           verbose = .true.)
  call Control_Mod_Starting_Time_Step_For_Swarm_Computation(first_dt_p,  &
                                                            verbose = .true.)

  ! Read physical models for each domain from control file
  do d = 1, n_dom
    call Control_Mod_Switch_To_Domain(d)  ! take proper control file
    call Read_Control_Physical_Models(flow(d), turb(d), mult(d), swarm(d))
  end do

  !----------------------------------------------------------!
  !   Allocate memory for all variables (over all domains)   !
  !----------------------------------------------------------!
  do d = 1, n_dom
    call Control_Mod_Switch_To_Domain(d)  ! take proper control file
    call Field_Mod_Allocate(flow(d), grid(d))
    call Turb_Mod_Allocate(turb(d), flow(d))
    call Swarm_Mod_Allocate(swarm(d), flow(d), turb(d))
    call Multiphase_Mod_Vof_Allocate(mult(d), flow(d))

    ! Read time step from root
    call Control_Mod_Switch_To_Root()
    call Control_Mod_Time_Step(flow(d) % dt, verbose=.true.)
    call Control_Mod_Switch_To_Domain(d)  ! go back to local domain's control

    ! Read numerical models from control file (after the memory is allocated)
    call Read_Control_Numerical(flow(d), turb(d), mult(d))

    call Grid_Mod_Find_Nodes_Cells(grid(d))
    call Grid_Mod_Calculate_Weights_Cells_To_Nodes(grid(d))  ! needed for front
    call Grid_Mod_Calculate_Global_Volumes(grid(d))
    call Field_Mod_Calculate_Grad_Matrix(flow(d))

    ! Allocate memory for linear systems of equations
    ! (You need face geomtry for this step)
    call Sol(d) % Create_Solver(grid(d))

    call Read_Control_Physical_Properties(flow(d), mult(d), swarm(d))
    call Read_Control_Boundary_Conditions(flow(d), turb(d), mult(d),   &
                                          turb_planes(d))
  end do

  ! Create interfaces
  call Control_Mod_Switch_To_Root()
  call Interface_Mod_Create(inter, grid, n_dom)

  ! First time step is one, unless read from backup otherwise
  first_dt = 0

  ! Read backup file if directed so, and set the "backup" to .true. or .false.
  do d = 1, n_dom
    call Control_Mod_Switch_To_Domain(d)  ! take proper control file
    call Backup_Mod_Load(flow(d), swarm(d), turb(d), mult(d),  &
                         time, first_dt, read_backup(d))

    ! Initialize variables
    if(.not. read_backup(d)) then
      call Initialize_Variables(flow(d), turb(d), mult(d), swarm(d), Sol(d))
      call Update_Boundary_Values(flow(d), turb(d), mult(d), 'ALL')
      call Results_Mod_Save_Surf(mult(d) % surf, curr_dt)
    end if

    if(mult(d) % model .eq. VOLUME_OF_FLUID) then
      if (read_backup(d))  then
        flow % piso_status = .false.
      end if
      call Multiphase_Mod_Vof_Physical_Properties(mult(d))
    end if

    ! Initialize monitoring points
    call Monitor(d) % Initialize(grid(d), read_backup(d), domain=d)

    ! Plane for calcution of overall mass fluxes
    call Control_Mod_Point_For_Monitoring_Planes(flow(d) % bulk % xp,  &
                                                 flow(d) % bulk % yp,  &
                                                 flow(d) % bulk % zp)

    ! Prepare ...
    call Bulk_Mod_Monitoring_Planes_Areas(flow(d) % bulk, grid(d))

    if( (turb(d) % model .eq. LES_SMAGORINSKY     .or.   &
         turb(d) % model .eq. HYBRID_LES_PRANDTL) .and.  &
         .not. read_backup(d)) then
      call Find_Nearest_Wall_Cell(turb(d))
    end if

    ! Print the areas of monitoring planes
    call Bulk_Mod_Print_Areas(flow(d) % bulk)

    ! Compute deltas for Spalart-Allmaras models
!   call Turb_Mod_Calculate_Deltas(turb(d))

  end do

  !---------------!
  !               !
  !   Time loop   !
  !               !
  !---------------!

  call Control_Mod_Switch_To_Root()
  call Control_Mod_Backup_Save_Interval  (backup % interval, verbose=.true.)
  call Control_Mod_Results_Save_Interval (result % interval, verbose=.true.)
  call Control_Mod_Save_Initial_Condition(result % initial,  verbose=.true.)
  if(mult(d) % model .eq. LAGRANGIAN_PARTICLES) then
    call Control_Mod_Swarm_Save_Interval(prsi, verbose=.true.)
  end if

  !-------------------------------------------------------------!
  !   Perform potential initialization in the first time step   !
  !-------------------------------------------------------------!
  if(first_dt .eq. 0) then
    do d = 1, n_dom
      call Control_Mod_Switch_To_Domain(d)  ! not sure if this call is needed
      call Control_Mod_Potential_Initialization(pot_init, .true.)
      if(pot_init) call Field_Mod_Potential_Initialization(flow(d), Sol(d))
    end do
  end if

  ! It will save results in .vtk or .cgns file format,
  ! depending on how the code was compiled
  ! First calls saves inside, second only the boundary cells
  if(result % initial) then
    do d = 1, n_dom
      call Results_Mod_Save(flow(d), turb(d), mult(d), swarm(d), first_dt,  &
                            plot_inside=.true., domain=d)
      call Results_Mod_Save(flow(d), turb(d), mult(d), swarm(d), first_dt,  &
                            plot_inside=.false., domain=d)
    end do
  end if

  do curr_dt = first_dt + 1, last_dt

    !------------------------------------!
    !   Preparations for new time step   !
    !------------------------------------!
    do d = 1, n_dom

      call Control_Mod_Switch_To_Domain(d)  ! not sure if this call is needed

      ! Update turbulent planes
      do tp = 1, turb_planes(d) % n_planes
        call Eddies_Mod_Superimpose(turb_planes(d) % plane(tp))
        call Eddies_Mod_Advance    (turb_planes(d) % plane(tp))
      end do

      if(d .eq. 1) time = time + flow(d) % dt

      ! Beginning of time step
      call User_Mod_Beginning_Of_Time_Step(flow(d), turb(d), mult(d),  &
                                           swarm(d), curr_dt, time)

      ! Start info boxes.
      call Info_Mod_Time_Start()
      call Info_Mod_Iter_Start()
      call Info_Mod_Bulk_Start()

      ! Initialize and print time info box
      if(d .eq. 1) then
        call Info_Mod_Time_Fill(curr_dt, time)
        call Info_Mod_Time_Print()
      end if

      ! Turbulence models initializations
      call Turb_Mod_Init(turb(d))

      ! Interface tracking
      if(mult(d) % model .eq. VOLUME_OF_FLUID) then
        call Multiphase_Mod_Vof_Main(mult(d), flow(d), turb(d), Sol(d), curr_dt)
        if(mult(d) % track_front) then
          call Results_Mod_Save_Front(mult(d) % Front, curr_dt)
!f_vs_s   call Results_Mod_Save_Surf (mult(d) % surf, curr_dt)
          call Results_Mod_Save(flow(d), turb(d), mult(d), swarm(d), curr_dt,  &
                                plot_inside=.true., domain=d)
        end if
        call Multiphase_Mod_Vof_Physical_Properties(mult(d))
      end if

      ! Lagrangian particle tracking
      if(mult(d) % model .eq. LAGRANGIAN_PARTICLES) then
        call User_Mod_Insert_Particles(flow(d), turb(d), mult(d),  &
                                       swarm(d), curr_dt, time)
      end if

    end do  ! through domains

    !--------------------------!
    !   Inner-iteration loop   !
    !--------------------------!
    call Control_Mod_Switch_To_Root()
    call Control_Mod_Max_Simple_Iterations(max_ini)
    call Control_Mod_Min_Simple_Iterations(min_ini)
    call Control_Mod_Tolerance_For_Simple_Algorithm(simple_tol)

    do ini = 1, max_ini

      ! Exchange data between domains
      call User_Mod_Interface_Exchange(inter, flow, n_dom)

      do d = 1, n_dom

        call Control_Mod_Switch_To_Domain(d)

        ! Beginning of iteration
        call User_Mod_Beginning_Of_Iteration(flow(d), turb(d), mult(d),  &
                                             swarm(d), curr_dt, time)

        call Info_Mod_Iter_Fill(ini)

        ! Compute velocity gradients
        call Field_Mod_Grad_Variable(flow(d), flow(d) % u)
        call Field_Mod_Grad_Variable(flow(d), flow(d) % v)
        call Field_Mod_Grad_Variable(flow(d), flow(d) % w)

        ! All three velocity components one after another
        call Compute_Momentum(flow(d), turb(d), mult(d), Sol(d), curr_dt, ini)
        call Compute_Pressure(flow(d), mult(d), Sol(d), curr_dt, ini)

        call Field_Mod_Calculate_Mass_Fluxes(flow(d), flow(d) % v_flux % n)
        call Correct_Velocity(flow(d), turb(d), mult(d), Sol(d), curr_dt, ini)

        call Piso_Algorithm(flow(d), turb(d), mult(d), Sol(d), ini)

        ! Energy (practically temperature)
        if(flow(d) % heat_transfer) then
          call Compute_Energy(flow(d), turb(d), mult(d), Sol(d), curr_dt, ini)
        end if

        ! Passive scalars
        do sc = 1, flow(d) % n_scalars
          call Compute_Scalar(flow(d), turb(d), mult(d), Sol(d),  &
                              curr_dt, ini, sc)
        end do

        ! Deal with turbulence (if you dare ;-))
        call Turb_Mod_Main(turb(d), Sol(d), curr_dt, ini)

        ! Update the values at boundaries
        call Convective_Outflow(flow(d), turb(d), mult(d))
        call Update_Boundary_Values(flow(d), turb(d), mult(d), 'ALL')

        ! End of the current iteration
        call Info_Mod_Iter_Print(d)

        ! End of iteration
        call User_Mod_End_Of_Iteration(flow(d), turb(d), mult(d), swarm(d),  &
                                       curr_dt, time)
      end do  ! through domains

      if(ini >= min_ini) then
        if( maxval(flow(1:n_dom) % u % res) <= simple_tol .and.  &
            maxval(flow(1:n_dom) % v % res) <= simple_tol .and.  &
            maxval(flow(1:n_dom) % w % res) <= simple_tol .and.  &
            maxval(flow(1:n_dom) % vol_res) <= simple_tol ) goto 1
      end if

    end do    ! through inner iterations

    !----------------------------------!
    !   End of the current time step   !
    !----------------------------------!
1   continue
    do d = 1, n_dom
      call Info_Mod_Bulk_Print(flow(d), d, n_dom)
    end do

    do d = 1, n_dom

      call Control_Mod_Switch_To_Domain(d)

      ! Write the values in monitoring points
      if(.not. flow(d) % heat_transfer) then
        call Monitor(d) % Write_4_Vars(curr_dt, flow(d))
      else
        call Monitor(d) % Write_5_Vars(curr_dt, flow(d))
      end if

      ! Calculate mean values
      call Turb_Mod_Calculate_Mean(turb(d), n_stat_t, curr_dt)
      call User_Mod_Calculate_Mean(turb(d), n_stat_t, curr_dt)

      ! Adjust pressure drops to keep the mass fluxes constant
      call Bulk_Mod_Adjust_P_Drops(flow(d) % bulk, flow(d) % dt)

      ! Lagrangian particle tracking
      if(mult(d) % model .eq. LAGRANGIAN_PARTICLES) then
        if(curr_dt >= first_dt_p) then
          call Swarm_Mod_Advance_Particles(swarm(d), curr_dt,  &
                                           n_stat_p, first_dt_p)
        end if
      end if

      ! Just before the end of time step
      call User_Mod_End_Of_Time_Step(flow(d), turb(d), mult(d), swarm(d),  &
                                     curr_dt, n_stat_t, n_stat_p, time)
    end do

    !----------------------!
    !   Save the results   !
    !----------------------!
    call Results_Mod_Main(curr_dt, last_dt, time, n_dom,  &
                          flow, turb, mult, swarm, exit_now)

    ! Ran more than a set wall clock time limit
    if(Info_Mod_Time_To_Exit() .or. exit_now) then
      goto 2
    end if

    ! Last time step reached; call user function for end of simulation
    if(curr_dt .eq. last_dt) then
      do d = 1, n_dom
        call Control_Mod_Switch_To_Domain(d)
        call User_Mod_End_Of_Simulation(flow(d), turb(d), mult(d), swarm(d),  &
                                        curr_dt, time)
      end do
    end if

  end do ! curr_dt until the last time step

2 if(this_proc < 2) print *, '# Exiting !'

  do d = 1, n_dom
    ! Close monitoring files
    call Monitor(d) % Finalize()

    ! Make the final call to user function
    call User_Mod_Before_Exit(grid(d))
  end do

  call Cpu_Timer % Stop('Main')
  call Cpu_Timer % Statistics

  !----------------------------!
  !   End parallel execution   !
  !----------------------------!
  call Comm_Mod_End

  end program

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
  use Interface_Mod
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
  character(len=7)      :: root_control    = 'control'
  character(len=9)      :: dom_control(MD) = 'control.n'
  integer               :: n, sc, tp
  real                  :: mass_res(MD)
  logical               :: backup(MD), save_now, exit_now
  type(Grid_Type)       :: grid(MD)        ! grid used in computations
  type(Field_Type)      :: flow(MD)        ! flow field we will be solving for
  type(Swarm_Type)      :: swarm(MD)       ! swarm of particles
  type(Surf_Type)       :: surf(MD)        ! interface between two phases
  type(Turb_Type)       :: turb(MD)        ! turbulence modelling
  type(Multiphase_Type) :: mult(MD)        ! multiphase modelling
  type(Solver_Type)     :: sol(MD)         ! linear solvers
  type(Turb_Plane_Type) :: turb_planes(MD) ! holder for synthetic turbulences
  type(Monitor_Type)    :: monitor(MD)     ! monitors
  type(Interface_Type)  :: inter(MD,MD)    ! interfaces between domains
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
  integer               :: n_dom           ! number of domains
  integer               :: d               ! domain counter
!==============================================================================!

  ! Initialize program profler
  call Cpu_Timer_Mod_Start('Main')

  ! Get starting time
  call system_clock(count_rate=sc_cr)  ! get system clock clock rate
  call system_clock(sc_ini)            ! get system clock initial count rate

  ! Initialize control file names
  root_control = 'control'             ! root control file name
  do d = 1, MD                         ! domain control file names
    write(dom_control(d), '(a8,i1)') 'control.', d
  end do

  ! Initialize variables
  time      =  0.      ! initialize time to zero
  backup(:) = .false.  ! can turn .true. in Load_Backup

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

  do d = 1, n_dom
    call Control_Mod_Switch_To_Domain(d)  ! take domain's d control file
    call Control_Mod_Problem_Name(problem_name(d))

    ! Load the finite volume grid
    call Grid_Mod_Load_Cns(grid(d), this_proc, domain=d)

    ! Allocate memory for communication
    call Grid_Mod_Allocate_Comm(grid(d))

    call Grid_Mod_Load_Geo(grid(d), this_proc, domain=d)
    call Grid_Mod_Create_Buffers(grid(d))
    call Grid_Mod_Create_Maps(grid(d))

    call Comm_Mod_Wait
  end do

  ! Create interfaces
  call Control_Mod_Switch_To_Root()
  call Interface_Mod_Create(inter, grid, n_dom)

  ! Allocate memory for working arrays
  call Work_Mod_Allocate(grid, rc=30, rf=1, rn=1, ic=4, if=0, in=1)

  ! Get the number of time steps from the control file
  call Control_Mod_Number_Of_Time_Steps(last_dt, verbose=.true.)
  call Control_Mod_Starting_Time_Step_For_Statistics(n_stat, verbose=.true.)

  ! Read physical models for each domain from control file
  do d = 1, n_dom
    call Control_Mod_Switch_To_Domain(d)  ! take proper control file
    call Read_Control_Physical(flow(d), turb(d), mult(d), swarm(d))
  end do

  ! Allocate memory for all variables (over all domains)
  do d = 1, n_dom
    call Control_Mod_Switch_To_Domain(d)  ! take proper control file
    call Field_Mod_Allocate(flow(d), grid(d))
    call Turb_Mod_Allocate(turb(d), flow(d))
    call Swarm_Mod_Allocate(swarm(d), flow(d))
    call Multiphase_Mod_Allocate(mult(d), flow(d))
    call User_Mod_Allocate(grid(d))

    ! Read numerical models from control file (after the memory is allocated)
    call Read_Control_Numerical(flow(d), turb(d), mult(d))

    call Grid_Mod_Calculate_Face_Geometry(grid(d))
    call Grid_Mod_Find_Nodes_Cells(grid(d))      ! for Lagrangian particle track
    call Grid_Mod_Find_Periodic_Faces(grid(d))   ! for Lagrangian particle track
    call Grid_Mod_Find_Cells_Faces(grid(d))      ! for Multiphase Module
    call Grid_Mod_Calculate_Global_Volumes(grid(d))

    ! Allocate memory for linear systems of equations
    ! (You need face geomtry for this step)
    call Solver_Mod_Create(sol(d), grid(d))

    call Load_Physical_Properties(flow(d), mult(d), swarm(d))

    call Load_Boundary_Conditions(flow(d), turb(d), mult(d), turb_planes(d))
  end do

  ! First time step is one, unless read from backup otherwise
  first_dt = 0

  ! Read backup file if directed so, and set the "backup" to .true. or .false.
  do d = 1, n_dom
    call Control_Mod_Switch_To_Domain(d)  ! take proper control file
    call Backup_Mod_Load(flow(d), swarm(d), turb(d), mult(d),  &
                         time, first_dt, n_stat, backup(d))

    ! Initialize variables
    if(.not. backup(d)) then
      call Initialize_Variables(flow(d), turb(d), mult(d), swarm(d))
      if(multiphase_model .eq. VOLUME_OF_FLUID) then
        call Multiphase_Mod_Update_Physical_Properties(mult(d))
      end if
    end if

    ! Initialize monitoring points
    call Monitor_Mod_Initialize(monitor(d), grid(d), backup(d), domain=d)

    ! Plane for calcution of overall mass fluxes
    call Control_Mod_Point_For_Monitoring_Planes(flow(d) % bulk % xp,  &
                                                 flow(d) % bulk % yp,  &
                                                 flow(d) % bulk % zp)

    ! Prepare ...
    call Bulk_Mod_Monitoring_Planes_Areas(flow(d) % bulk, grid(d))

    if(turb(d) % model .eq. LES_SMAGORINSKY .and. .not. backup(d)) then
      call Find_Nearest_Wall_Cell(turb(d))
    end if

    if(turb(d) % model .eq. HYBRID_LES_PRANDTL .and. .not. backup(d)) then
      call Find_Nearest_Wall_Cell(turb(d))
    end if

    ! Prepare the gradient matrix for velocities
    call Field_Mod_Calculate_Grad_Matrix(flow(d))

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
  call Control_Mod_Backup_Save_Interval (bsi,    verbose=.true.)
  call Control_Mod_Results_Save_Interval(rsi,    verbose=.true.)
  call Control_Mod_Wall_Time_Max_Hours  (wt_max, verbose=.true.)
  wt_max = wt_max * 3600  ! make it in seconds

  ! It will save results in .vtk or .cgns file format,
  ! depending on how the code was compiled
  ! First calls saves inside, second only the boundary cells
  do d = 1, n_dom
    call Save_Results(flow(d), turb(d), mult(d), swarm(d), first_dt,  &
                      plot_inside=.true., domain=d)
    call Save_Results(flow(d), turb(d), mult(d), swarm(d), first_dt,  &
                      plot_inside=.false., domain=d)
    call Save_Swarm(swarm(d), first_dt, domain=d)
  end do

  do n = first_dt + 1, last_dt

    !------------------------------------!
    !   Preparations for new time step   !
    !------------------------------------!
    do d = 1, n_dom

      call Control_Mod_Switch_To_Root()  ! read time step from root
      call Control_Mod_Time_Step(flow(d) % dt, verbose=.true.)

      call Control_Mod_Switch_To_Domain(d)  ! not sure if this call is needed

      ! Update turbulent planes
      do tp = 1, turb_planes(d) % n_planes
        call Eddies_Mod_Superimpose(turb_planes(d) % plane(tp))
        call Eddies_Mod_Advance    (turb_planes(d) % plane(tp))
      end do

      if(d .eq. 1) time = time + flow(d) % dt

      ! Beginning of time step
      call User_Mod_Beginning_Of_Time_Step(flow(d), turb(d), mult(d),  &
                                           swarm(d), n, time)

      ! Start info boxes.
      call Info_Mod_Time_Start()
      call Info_Mod_Iter_Start()
      call Info_Mod_Bulk_Start()

      ! Initialize and print time info box
      call system_clock(sc_cur)
      if(d .eq. 1) then
        call Info_Mod_Time_Fill( n, time, real(sc_cur-sc_ini)/real(sc_cr) )
        call Info_Mod_Time_Print()
      end if

      ! Turbulence models initializations
      call Turb_Mod_Init(turb(d), sol(d))

      ! Volume of Fluid
      if(multiphase_model .eq. VOLUME_OF_FLUID) then
        flow(d) % m_flux % o = flow(d) % m_flux % n / flow(d) % density_f
        ! Update the values at boundaries
        call Update_Boundary_Values(flow(d), turb(d), mult(d))
        call Multiphase_Mod_Compute_Vof(mult(d), sol(d), flow(d) % dt, n)
      else
        flow(d) % m_flux % o = flow(d) % m_flux % n
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
      call Interface_Mod_Exchange(inter, flow, n_dom)

      do d = 1, n_dom

        call Control_Mod_Switch_To_Domain(d)

        ! Beginning of iteration
        call User_Mod_Beginning_Of_Iteration(flow(d), turb(d), mult(d),  &
                                             swarm(d), n, time)

        call Info_Mod_Iter_Fill(ini)

        call Field_Mod_Grad_Pressure(flow(d), flow(d) % p,  &
                                     flow(d) % density,     &
                                     grav_x, grav_y, grav_z)

        ! Compute velocity gradients
        call Field_Mod_Grad_Variable(flow(d), flow(d) % u)
        call Field_Mod_Grad_Variable(flow(d), flow(d) % v)
        call Field_Mod_Grad_Variable(flow(d), flow(d) % w)

        ! All velocity components one after another
        call Compute_Momentum(flow(d), turb(d), mult(d), 1, sol(d),  &
                              flow(d) % dt, ini)
        call Compute_Momentum(flow(d), turb(d), mult(d), 2, sol(d),  &
                              flow(d) % dt, ini)
        call Compute_Momentum(flow(d), turb(d), mult(d), 3, sol(d),  &
                              flow(d) % dt, ini)

        ! Refresh buffers for a % sav before discretizing for pressure
        ! (Can this call be somewhere in Compute Pressure?)
        call Grid_Mod_Exchange_Real(grid(d), sol(d) % a % sav)

        call Balance_Mass(flow(d))
        call Compute_Pressure(flow(d), mult(d), sol(d), flow(d) % dt, ini)

        call Field_Mod_Grad_Pressure_Correction(flow(d), flow(d) % pp)

        call Field_Mod_Calculate_Fluxes(flow(d), flow(d) % m_flux % n)
        mass_res(d) = Correct_Velocity(flow(d), sol(d), flow(d) % dt, ini)

        ! Energy (practically temperature)
        if(heat_transfer) then
          call Compute_Energy(flow(d), turb(d), mult(d), sol(d),  &
                              flow(d) % dt, ini)
        end if

        ! Passive scalars
        do sc = 1, flow(d) % n_scalars
          call Compute_Scalar(flow(d), turb(d), mult(d), sol(d),  &
                              flow(d) % dt, ini, sc)
        end do

        ! Deal with turbulence (if you dare ;-))
        call Turb_Mod_Main(turb(d), sol(d), n, ini)

        ! Update the values at boundaries
        call Convective_Outflow(flow(d), turb(d), mult(d), flow(d) % dt)
        call Update_Boundary_Values(flow(d), turb(d), mult(d))

        ! End of the current iteration
        call Info_Mod_Iter_Print(d)

        ! End of iteration
        call User_Mod_End_Of_Iteration(flow(d), turb(d), mult(d), swarm(d),  &
                                       n, time)
      end do  ! through domains

      if(ini >= min_ini) then
        if( maxval(flow(1:n_dom) % u  % res) <= simple_tol .and.  &
            maxval(flow(1:n_dom) % v  % res) <= simple_tol .and.  &
            maxval(flow(1:n_dom) % w  % res) <= simple_tol .and.  &
            maxval(mass_res(1:n_dom))        <= simple_tol ) goto 1
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
      if(.not. heat_transfer) then
        call Monitor_Mod_Write_4_Vars(monitor(d), n, flow(d))
      else
        call Monitor_Mod_Write_5_Vars(monitor(d), n, flow(d))
      end if

      ! Calculate mean values
      call Turb_Mod_Calculate_Mean(turb(d), n_stat, n)
      call User_Mod_Calculate_Mean(turb(d), n_stat, n)

      ! Adjust pressure drops to keep the mass fluxes constant
      call Bulk_Mod_Adjust_P_Drops(flow(d) % bulk, flow(d) % dt)

      ! Just before the end of time step
      call User_Mod_End_Of_Time_Step(flow(d), turb(d), mult(d), swarm(d),  &
                                     n, time)
    end do

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
      do d = 1, n_dom
        call Control_Mod_Switch_To_Domain(d)
        call Backup_Mod_Save(flow(d), swarm(d), turb(d), mult(d),  &
                             time, n, n_stat, domain=d)
      end do
    end if

    ! Is it time to save results for post-processing
    if(save_now           .or.  &
       exit_now           .or.  &
       mod(n, rsi) .eq. 0 .or.  &
       real(sc_cur-sc_ini)/real(sc_cr) > wt_max) then

      do d = 1, n_dom
        call Control_Mod_Switch_To_Domain(d)
        call Save_Results(flow(d), turb(d), mult(d), swarm(d), n,  &
                          plot_inside=.true., domain=d)
        call Save_Results(flow(d), turb(d), mult(d), swarm(d), n,  &
                          plot_inside=.false., domain=d)
        call Save_Swarm(swarm(d), n)

        if(multiphase_model .eq. VOLUME_OF_FLUID) then
!         call Surf_Mod_Allocate(surf(d), flow(d))
!         call Surf_Mod_Place_At_Var_Value(surf(d),        &
!                                          mult(d) % vof,  &
!                                          sol(d),         &
!                                          0.5,            &
!                                          .false.)  ! don't print messages
!         call Save_Surf(surf(d), n)
!         call Surf_Mod_Clean(surf(d))
        end if

        ! Write results in user-customized format
        call User_Mod_Save_Results(flow(d), turb(d), mult(d), swarm(d), n)
        ! call User_Mod_Save_Swarm(swarm(d), n)  to be done!

      end do  ! through domains
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
  do d = 1, n_dom
    call Control_Mod_Switch_To_Domain(d)

    call User_Mod_End_Of_Simulation(flow(d), turb(d), mult(d), swarm(d),  &
                                    n, time)

    ! Save backup and post-processing files at exit
    call Backup_Mod_Save(flow(d), swarm(d), turb(d), mult(d),  &
                         time, n, n_stat, domain=d)
    call Save_Results(flow(d), turb(d), mult(d), swarm(d), n,  &
                      plot_inside=.true., domain=d)
    call Save_Results(flow(d), turb(d), mult(d), swarm(d), n,  &
                      plot_inside=.false., domain=d)

    ! Write results in user-customized format
    call User_Mod_Save_Results(flow(d), turb(d), mult(d), swarm(d), n)
  end do

  if(this_proc < 2) then
    open(9, file='stop')
    close(9)
  end if

2 if(this_proc  < 2) print *, '# Exiting !'

  do d = 1, n_dom
    ! Close monitoring files
    call Monitor_Mod_Finalize(monitor(d))

    ! Make the final call to user function
    call User_Mod_Before_Exit(grid(d))
  end do

  call Cpu_Timer_Mod_Stop('Main')
  call Cpu_Timer_Mod_Statistics

  !----------------------------!
  !   End parallel execution   !
  !----------------------------!
  call Comm_Mod_End

  end program

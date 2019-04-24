!==============================================================================!
  program Processor
!------------------------------------------------------------------------------!
!   Unstructured finite volume 'LES'/RANS solver.                              !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Const_Mod
  use Comm_Mod
  use Name_Mod,    only: problem_name
  use Field_Mod,   only: Field_Type, Field_Mod_Allocate, heat_transfer
  use Turb_Mod
  use Grid_Mod
  use Grad_Mod
  use Bulk_Mod
  use Var_Mod,     only: Var_Type
  use Solver_Mod,  only: Solver_Mod_Create
  use Info_Mod
  use Work_Mod,    only: Work_Mod_Allocate
  use User_Mod
  use Save_Results_Mod
  use Control_Mod
  use Monitor_Mod
  use Backup_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Calling]------------------------------------!
  real :: Correct_Velocity
!----------------------------------[Locals]------------------------------------!
  integer           :: n, sc
  real              :: mass_res, wall_time_start, wall_time_current
  character(len=80) :: name_save
  logical           :: backup, save_now, exit_now
  type(Grid_Type)   :: grid        ! grid used in computations
  type(Field_Type)  :: flow        ! flow field we will be solving for
  type(Swarm_Type)  :: swarm       ! swarm of particles
  type(Turb_Type)   :: turb        ! turbulence modelling
  type(Solver_Type) :: sol         ! linear solvers
  real              :: time        ! physical time of the simulation
  integer           :: first_dt    ! first time step in this run
  integer           :: last_dt     ! number of time steps
  integer           :: max_ini     ! max number of inner iterations
  integer           :: min_ini     ! min number of inner iterations
  integer           :: n_stat      ! starting time step for statistic
  integer           :: ini         ! inner iteration counter
  integer           :: bsi, rsi    ! backup and results save interval
  real              :: simple_tol  ! tolerance for SIMPLE algorithm
!==============================================================================!

  ! Get starting time
  call cpu_time(wall_time_start)
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
  call Control_Mod_Open_File()
  call Control_Mod_Problem_Name(problem_name)

  ! Load the finite volume grid
  call Grid_Mod_Load_Cns(grid, this_proc)

  ! Allocate memory for working arrays and comm.
  call Work_Mod_Allocate(grid, 30, 1, 1, 4)
  call Comm_Mod_Allocate(grid)

  call Grid_Mod_Load_Geo(grid, this_proc)
  call Comm_Mod_Create_Buffers(grid)
  call Comm_Mod_Load_Maps(grid)       ! maps should move to .cns file soon

  call Comm_Mod_Wait

  ! Get the number of time steps from the control file
  call Control_Mod_Number_Of_Time_Steps(last_dt, verbose=.true.)
  call Control_Mod_Starting_Time_Step_For_Statistics(n_stat, verbose=.true.)

  ! Read physical models from control file
  call Read_Control_Physical(flow, swarm, turb, backup)

  ! Read numerical models from control file
  call Read_Control_Numerical(flow, turb)

  ! Allocate memory for all variables
  call Field_Mod_Allocate(flow, grid)
  call Grad_Mod_Allocate(grid)
  call Turb_Mod_Allocate(turb, flow)
  call Swarm_Mod_Allocate(swarm, flow)
  call User_Mod_Allocate(grid)

  call Grid_Mod_Calculate_Face_Geometry(grid)
  call Grid_Mod_Find_Nodes_Cells(grid)         ! for Lagrangian particle track

  ! Allocate memory for linear systems of equations
  ! (You need face geomtry for this step)
  call Solver_Mod_Create(sol, grid)

  call Load_Physical_Properties(grid)

  call Load_Boundary_Conditions(flow, turb, backup)

  ! First time step is one, unless read from backup otherwise
  first_dt = 0

  ! Read backup file if directed so, and set the "backup" to .true. or .false.
  call Backup_Mod_Load(flow, swarm, turb, first_dt, n_stat, backup) 

  ! Initialize variables
  if(.not. backup) then
    call Initialize_Variables(flow, turb)
    call Comm_Mod_Wait
  end if

  ! Initialize monitoring points
  call Monitor_Mod_Initialize(grid, backup)

  ! Plane for calcution of overall mass fluxes
  call Control_Mod_Point_For_Monitoring_Planes(flow % bulk % xp,  &
                                               flow % bulk % yp,  &
                                               flow % bulk % zp)

  ! Prepare ...
  call Bulk_Mod_Monitoring_Planes_Areas(flow % bulk, grid)
  call Grad_Mod_Find_Bad_Cells         (grid)

  if(turbulence_model .eq. LES_SMAGORINSKY .and. .not. backup) then
    call Find_Nearest_Wall_Cell(grid)
  end if

  ! Prepare the gradient matrix for velocities
  call Compute_Gradient_Matrix(grid, .true.)

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
  call Control_Mod_Backup_Save_Interval(bsi, verbose=.true.)
  call Control_Mod_Results_Save_Interval(rsi, verbose=.true.)

  ! Form the file name before the time, in case ...
  ! ... user just wants to save files from backup
  name_save = problem_name
  write(name_save(len_trim(problem_name)+1:                    &
                  len_trim(problem_name)+3), '(a3)')   '-ts'
  write(name_save(len_trim(problem_name)+4:                    &
                  len_trim(problem_name)+9), '(i6.6)') first_dt

  ! It will save results in .vtk or .cgns file format,
  ! depending on how the code was compiled
  call Save_Results(flow, turb, problem_name)
  call Save_Swarm (swarm, problem_name)

  do n = first_dt + 1, last_dt

    time = time + flow % dt

    ! Beginning of time steo
    call User_Mod_Beginning_Of_Time_Step(flow, turb, swarm, n, time)

    ! Start info boxes.
    call Info_Mod_Time_Start()
    call Info_Mod_Iter_Start()
    call Info_Mod_Bulk_Start()

    ! Initialize and print time info box
    call cpu_time(wall_time_current)
    call Info_Mod_Time_Fill( n, time, (wall_time_current-wall_time_start) )
    call Info_Mod_Time_Print()

    if(turbulence_model .eq. DES_SPALART) then
      call Calculate_Shear_And_Vorticity(flow)
      call Calculate_Vorticity (flow)
    end if

    if(turbulence_model .eq. LES_SMAGORINSKY .or.  &
       turbulence_model .eq. LES_DYNAMIC     .or.  &
       turbulence_model .eq. LES_WALE) then
      call Calculate_Shear_And_Vorticity(flow)
      if(turbulence_model .eq. LES_DYNAMIC) then
        call Turb_Mod_Vis_T_Dynamic(turb, sol)
      end if
      if(turbulence_model .eq. LES_WALE) then
        call Turb_Mod_Vis_T_Wale(turb)
      end if
      call Turb_Mod_Vis_T_Smagorinsky(turb)
    end if

    if(turbulence_model .eq. HYBRID_LES_RANS) then
      call Turb_Mod_Vis_T_Dynamic(turb, sol)
      call Turb_Mod_Vis_T_Hybrid (turb)
    end if

    call Convective_Outflow(flow, flow % dt)
    if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
       turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
      call Turb_Mod_Vis_T_Rsm(turb)
    end if

    !--------------------------!
    !   Inner-iteration loop   !
    !--------------------------!
    call Control_Mod_Max_Simple_Iterations(max_ini)
    call Control_Mod_Min_Simple_Iterations(min_ini)

    do ini = 1, max_ini

      call Info_Mod_Iter_Fill(ini)

      call Grad_Mod_Pressure(flow % p)

      ! Compute velocity gradients
      call Grad_Mod_Variable(flow % u, .true.)
      call Grad_Mod_Variable(flow % v, .true.)
      call Grad_Mod_Variable(flow % w, .true.)

      ! All velocity components one after another
      call Compute_Momentum(flow, turb, 1, sol, flow % dt, ini)
      call Compute_Momentum(flow, turb, 2, sol, flow % dt, ini)
      call Compute_Momentum(flow, turb, 3, sol, flow % dt, ini)

      ! Refresh buffers for a % sav before discretizing for pressure
      ! (Can this call be somewhere in Compute Pressure?)
      call Comm_Mod_Exchange_Real(grid, sol % a % sav)

      call Balance_Mass(flow)
      call Compute_Pressure(flow, sol, flow % dt, ini)

      call Grad_Mod_Pressure(flow % pp)

      call Bulk_Mod_Calculate_Fluxes(grid, flow % bulk, flow % flux)
      mass_res = Correct_Velocity(flow, sol, flow % dt, ini)

      ! Energy (practically temperature)
      if(heat_transfer) then
        call Compute_Energy(flow, turb, sol, flow % dt, ini)
      end if

      ! Passive scalars
      do sc = 1, flow % n_scalars
        call Compute_Scalar(flow, turb, sol, flow % dt, ini, sc)
      end do

      if(turbulence_model .eq. K_EPS) then

        ! Update the values at boundaries
        call Update_Boundary_Values(flow, turb)

        call Calculate_Shear_And_Vorticity(flow)

        call Turb_Mod_Compute_Variable(turb, sol, ini, turb % kin, n)
        call Turb_Mod_Compute_Variable(turb, sol, ini, turb % eps, n)

        call Turb_Mod_Vis_T_K_Eps(turb)

      end if

      if(turbulence_model .eq. K_EPS_ZETA_F .or. &
         turbulence_model .eq. HYBRID_LES_RANS) then
        call Calculate_Shear_And_Vorticity(flow)

        call Turb_Mod_Compute_Variable(turb, sol, ini, turb % kin, n)
        call Turb_Mod_Compute_Variable(turb, sol, ini, turb % eps, n)

        if(heat_transfer) then
          call Turb_Mod_Calculate_Heat_Flux(turb)
          call Turb_Mod_Compute_Variable(turb, sol, ini, turb % t2, n)
        end if

        call Update_Boundary_Values(flow, turb)

        call Turb_Mod_Compute_F22(turb, sol, ini, turb % f22)
        call Turb_Mod_Compute_Variable(turb, sol, ini, turb % zeta, n)

        call Turb_Mod_Vis_T_K_Eps_Zeta_F(turb)

      end if

      if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
         turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then

        ! Update the values at boundaries
        call Update_Boundary_Values(flow, turb)

        call Time_And_Length_Scale(grid, turb)

        call Grad_Mod_Variable(flow % u, .true.)
        call Grad_Mod_Variable(flow % v, .true.)
        call Grad_Mod_Variable(flow % w, .true.)

        call Turb_Mod_Compute_Stress(turb, sol, ini, turb % uu, n)
        call Turb_Mod_Compute_Stress(turb, sol, ini, turb % vv, n)
        call Turb_Mod_Compute_Stress(turb, sol, ini, turb % ww, n)

        call Turb_Mod_Compute_Stress(turb, sol, ini, turb % uv, n)
        call Turb_Mod_Compute_Stress(turb, sol, ini, turb % uw, n)
        call Turb_Mod_Compute_Stress(turb, sol, ini, turb % vw, n)

        if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
          call Turb_Mod_Compute_F22(turb, sol, ini, turb % f22)
        end if

        call Turb_Mod_Compute_Stress(turb, sol, ini, turb % eps, n)

        call Turb_Mod_Vis_T_Rsm(turb)

        if(heat_transfer) then
          call Turb_Mod_Calculate_Heat_Flux(turb)
        end if
      end if

      if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
         turbulence_model .eq. DES_SPALART) then
        call Calculate_Shear_And_Vorticity(flow)
        call Calculate_Vorticity(flow)

        ! Update the values at boundaries
        call Update_Boundary_Values(flow, turb)

        call Turb_Mod_Compute_Variable(turb, sol, ini, turb % vis, n)
        call Turb_Mod_Vis_T_Spalart_Allmaras(turb)
      end if

      ! Update the values at boundaries
      call Update_Boundary_Values(flow, turb)

      ! End of the current iteration
      call Info_Mod_Iter_Print()

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
    call User_Mod_Calculate_Mean(flow, n_stat, n)

    ! Adjust pressure drops to keep the mass fluxes constant
    call Bulk_Mod_Adjust_P_Drops(flow % bulk, flow % dt)

    !----------------------!
    !   Save the results   !
    !----------------------!
    inquire(file='exit_now', exist=exit_now)
    inquire(file='save_now', exist=save_now)

    ! Form the file name
    name_save = problem_name
    write(name_save(len_trim(problem_name)+1:                    &
                    len_trim(problem_name)+3), '(a3)')   '-ts'
    write(name_save(len_trim(problem_name)+4:                    &
                    len_trim(problem_name)+9), '(i6.6)') n

    ! Is it time to save the backup file?
    if(save_now .or. exit_now .or. mod(n, bsi) .eq. 0) then
      call Backup_Mod_Save(flow, swarm, turb, n, n_stat, name_save)
    end if

    ! Is it time to save results for post-processing
    if(save_now .or. exit_now .or. mod(n, rsi) .eq. 0) then
      call Comm_Mod_Wait
      call Save_Results(flow, turb, name_save)
      call Save_Swarm (swarm, name_save)

      ! Write results in user-customized format
      call User_Mod_Save_Results(flow, turb, name_save)
      ! call User_Mod_Save_Swarm(swarm, name_save)  to be done!
    end if

    ! Just before the end of time step
    call User_Mod_End_Of_Time_Step(flow, turb, swarm, n, time)

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

  end do ! n, number of time steps

  ! After the time loop, decrease "n" since ...
  ! ... it is one above the loop boundaries here
  n = n - 1

  ! Save backup and post-processing files at exit
  call Comm_Mod_Wait
  call Save_Results(flow, turb, name_save)
  call Save_Swarm (swarm, name_save)
  call Backup_Mod_Save(flow, swarm, turb, n, n_stat, name_save)

  ! Write results in user-customized format
  call User_Mod_Save_Results(flow, turb, name_save)
  ! call User_Mod_Save_Swarm(swarm, name_save)  to be done!

  if(this_proc < 2) then
    open(9, file='stop')
    close(9)
  end if

2 if(this_proc  < 2) print *, '# Exiting !'

  ! Close monitoring files
  call Monitor_Mod_Finalize

  ! Make the final call to user function
  call User_Mod_Before_Exit(grid)

  !----------------------------!
  !   End parallel execution   !
  !----------------------------!
  call Comm_Mod_End

  end program

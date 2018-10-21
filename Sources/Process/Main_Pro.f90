!==============================================================================!
  program Processor
!------------------------------------------------------------------------------!
!   Unstructured finite volume 'LES'/RANS solver.                              !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Name_Mod, only: problem_name
  use Const_Mod
  use Flow_Mod
  use Les_Mod
  use Comm_Mod
  use Rans_Mod
  use Grid_Mod
  use Grad_Mod
  use Bulk_Mod
  use Var_Mod
  use Solvers_Mod, only: d
  use Info_Mod
  use User_Mod
  use Control_Mod
  use Monitor_Mod
  use Backup_Mod
  use Numerics_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Calling]------------------------------------!
  real :: Correct_Velocity
!----------------------------------[Locals]------------------------------------!
  integer           :: n, us
  real              :: mass_res, wall_time_start, wall_time_current
  character(len=80) :: name_save
  logical           :: backup, save_now, exit_now
  type(Grid_Type)   :: grid        ! grid used in computations
  real              :: time        ! physical time
  real              :: dt          ! time step
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
  call Load_Cns(grid, this_proc)

  call Allocate_Memory(grid)
  call Load_Geo(grid, this_proc)
  call Comm_Mod_Create_Buffers(grid)
  call Comm_Mod_Load_Maps(grid)     ! Maps should move to .cns file soon

  ! This is actually pretty bad - this command should be in Load_Geo
  call Comm_Mod_Exchange_Real(grid, grid % vol(-grid % n_bnd_cells))

  call Matrix_Mod_Topology(grid, a)
  call Matrix_Mod_Topology(grid, d)

  call Comm_Mod_Wait

  ! Get the number of time steps from the control file
  call Control_Mod_Number_Of_Time_Steps(last_dt, verbose=.true.)
  call Control_Mod_Starting_Time_Step_For_Statistics(n_stat, verbose=.true.)

  call Allocate_Variables(grid)

  call Calculate_Face_Geometry(grid)
  call Load_Physical_Properties(grid)

  call Load_Boundary_Conditions(grid, backup)

  ! First time step is one, unless read from backup otherwise
  first_dt = 0
  call Backup_Mod_Load(grid, first_dt, n_stat, backup)  ! "backup" can change

  ! Read physical models from control file
  call Read_Physical(grid, backup)

  ! Initialize variables
  if(.not. backup) then
    call Initialize_Variables(grid)
    call Comm_Mod_Wait
  end if

  ! Retreive information about monitoring points from control file
  call Control_Mod_Monitoring_Points(.true.)

  ! Initialize monitoring points
  call Monitor_Mod_Initialize(grid, backup)

  ! Plane for calcution of overall mass fluxes
  call Control_Mod_Point_For_Monitoring_Planes(bulk % xp,  &
                                               bulk % yp,  &
                                               bulk % zp)

  ! Prepare ...
  call Bulk_Mod_Monitoring_Planes_Areas(grid, bulk)
  call Grad_Mod_Find_Bad_Cells         (grid)

  if(turbulence_model .eq. LES_SMAGORINSKY .and. .not. backup) then
    call Find_Nearest_Wall_Cell(grid)
  end if

  ! Prepare the gradient matrix for velocities
  call Compute_Gradient_Matrix(grid, .true.)

  ! Print the areas of monitoring planes
  if(this_proc < 2) then
    write(*,'(a7,es12.5)') ' # Ax :', bulk % area_x
    write(*,'(a7,es12.5)') ' # Ay :', bulk % area_y
    write(*,'(a7,es12.5)') ' # Az :', bulk % area_z
  end if

  !---------------!
  !               !
  !   Time loop   !
  !               !
  !---------------!

  call Control_Mod_Time_Step(dt, verbose=.true.)
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
  call Save_Results(grid, problem_name)

  do n = first_dt + 1, last_dt

    time = time + dt

    ! Beginning of time steo
    call User_Mod_Beginning_Of_Time_Step(grid, n, time)

    ! Start info boxes.
    call Info_Mod_Time_Start()
    call Info_Mod_Iter_Start()
    call Info_Mod_Bulk_Start()

    ! Initialize and print time info box
    call cpu_time(wall_time_current)
    call Info_Mod_Time_Fill( n, time, (wall_time_current-wall_time_start) )
    call Info_Mod_Time_Print()

    if(turbulence_model .eq. DES_SPALART) then
      call Calculate_Shear_And_Vorticity(grid)
      call Calculate_Vorticity (grid)
    end if

    if(turbulence_model .eq. LES_SMAGORINSKY .or.  &
       turbulence_model .eq. LES_DYNAMIC     .or.  &
       turbulence_model .eq. LES_WALE) then
      call Calculate_Shear_And_Vorticity(grid)
      if(turbulence_model .eq. LES_DYNAMIC) call Calculate_Sgs_Dynamic(grid)
      if(turbulence_model .eq. LES_WALE)    call Calculate_Sgs_Wale(grid)
      call Calculate_Sgs(grid)
    end if

    If(turbulence_model .eq. HYBRID_LES_RANS) then
      call Calculate_Sgs_Dynamic(grid)
      call Calculate_Sgs_Hybrid(grid)
    end if

    call Convective_Outflow(grid, dt)
    if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
       turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
      call Calculate_Vis_T_Rsm(grid)
    end if

    !--------------------------!
    !   Inner-iteration loop   !
    !--------------------------!
    call Control_Mod_Max_Simple_Iterations(max_ini)
    call Control_Mod_Min_Simple_Iterations(min_ini)

    do ini = 1, max_ini

      call Info_Mod_Iter_Fill(ini)

      call Grad_Mod_For_P(grid, p % n, p % x, p % y, p % z)

      ! Compute velocity gradients
      call Grad_Mod_For_Phi(grid, u % n, 1, u % x, .true.)
      call Grad_Mod_For_Phi(grid, u % n, 2, u % y, .true.)
      call Grad_Mod_For_Phi(grid, u % n, 3, u % z, .true.)
      call Grad_Mod_For_Phi(grid, v % n, 1, v % x, .true.)
      call Grad_Mod_For_Phi(grid, v % n, 2, v % y, .true.)
      call Grad_Mod_For_Phi(grid, v % n, 3, v % z, .true.)
      call Grad_Mod_For_Phi(grid, w % n, 1, w % x, .true.)
      call Grad_Mod_For_Phi(grid, w % n, 2, w % y, .true.)
      call Grad_Mod_For_Phi(grid, w % n, 3, w % z, .true.)

      ! u velocity component
      call Compute_Momentum(grid, dt, ini, u,      &
                  u % x,     u % y,     u % z,     &
                  grid % sx, grid % sy, grid % sz, &
                  grid % dx, grid % dy, grid % dz, &
                  p % x,     v % x,     w % x)     ! dp/dx, dv/dx, dw/dx

      ! v velocity component
      call Compute_Momentum(grid, dt, ini, v,      &
                  v % y,     v % x,     v % z,     &
                  grid % sy, grid % sx, grid % sz, &
                  grid % dy, grid % dx, grid % dz, &
                  p % y,     u % y,     w % y)     ! dp/dy, du/dy, dw/dy

      ! w velocity component
      call Compute_Momentum(grid, dt, ini, w,      &
                  w % z,     w % x,     w % y,     &
                  grid % sz, grid % sx, grid % sy, &
                  grid % dz, grid % dx, grid % dy, &
                  p % z,     u % z,     v % z)     ! dp/dz, du/dz, dv/dz

      ! Refresh buffers for a % sav before discretizing for pressure
      call Comm_Mod_Exchange_Real(grid, a % sav)

      call Balance_Mass(grid)
      call Compute_Pressure_Simple(grid, dt, ini)

      call Grad_Mod_For_P(grid,  pp % n, p % x, p % y, p % z)

      call Bulk_Mod_Compute_Fluxes(grid, bulk, flux)
      mass_res = Correct_Velocity(grid, dt, ini) !  project the velocities

      ! Energy (practically temperature)
      if(heat_transfer) then
        call Compute_Energy(grid, dt, ini, t)
      end if

      ! User scalars
      do us = 1, n_user_scalars
        call User_Mod_Compute_Scalar(grid, dt, ini, user_scalar(us))
      end do

      ! Rans models
      if(turbulence_model .eq. K_EPS) then

        ! Update the values at boundaries
        call Update_Boundary_Values(grid)

        call Calculate_Shear_And_Vorticity(grid)

        call Compute_Turbulent(grid, dt, ini, kin, n)
        call Compute_Turbulent(grid, dt, ini, eps, n)

        call Calculate_Vis_T_K_Eps(grid)

        if(heat_transfer) then
          call Calculate_Heat_Flux(grid)
        end if
      end if

      if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
         turbulence_model .eq. HYBRID_LES_RANS) then
        call Calculate_Shear_And_Vorticity(grid)

        call Compute_Turbulent(grid, dt, ini, kin, n)
        call Compute_Turbulent(grid, dt, ini, eps, n)
        call Update_Boundary_Values(grid)

        call Compute_F22(grid, ini, f22)
        call Compute_Turbulent(grid, dt, ini, zeta, n)

        call Calculate_Vis_T_K_Eps_Zeta_F(grid)

        if(heat_transfer) then
          call Calculate_Heat_Flux(grid)
        end if
      end if

      if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
         turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then

        ! Update the values at boundaries
        call Update_Boundary_Values(grid)

        if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
          call Time_And_Length_Scale(grid)
        end if

        call Grad_Mod_For_Phi(grid, u % n, 1, u % x,.true.)    ! dU/dx
        call Grad_Mod_For_Phi(grid, u % n, 2, u % y,.true.)    ! dU/dy
        call Grad_Mod_For_Phi(grid, u % n, 3, u % z,.true.)    ! dU/dz

        call Grad_Mod_For_Phi(grid, v % n, 1, v % x,.true.)    ! dV/dx
        call Grad_Mod_For_Phi(grid, v % n, 2, v % y,.true.)    ! dV/dy
        call Grad_Mod_For_Phi(grid, v % n, 3, v % z,.true.)    ! dV/dz

        call Grad_Mod_For_Phi(grid, w % n, 1, w % x,.true.)    ! dW/dx
        call Grad_Mod_For_Phi(grid, w % n, 2, w % y,.true.)    ! dW/dy
        call Grad_Mod_For_Phi(grid, w % n, 3, w % z,.true.)    ! dW/dz

        call Compute_Stresses(grid, dt, ini, uu, n)
        call Compute_Stresses(grid, dt, ini, vv, n)
        call Compute_Stresses(grid, dt, ini, ww, n)

        call Compute_Stresses(grid, dt, ini, uv, n)
        call Compute_Stresses(grid, dt, ini, uw, n)
        call Compute_Stresses(grid, dt, ini, vw, n)

        if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
          call Compute_F22(grid, ini, f22)
        end if

        call Compute_Stresses(grid, dt, ini, eps, n)

        call Calculate_Vis_T_Rsm(grid)

        if(heat_transfer) then
          call Calculate_Heat_Flux(grid)
        end if
      end if

      if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
         turbulence_model .eq. DES_SPALART) then
        call Calculate_Shear_And_Vorticity(grid)
        call Calculate_Vorticity(grid)

        ! Update the values at boundaries
        call Update_Boundary_Values(grid)

        call Compute_Turbulent(grid, dt, ini, vis, n)
        call Calculate_Vis_T_Spalart_Allmaras(grid)
      end if

      ! Update the values at boundaries
      call Update_Boundary_Values(grid)

      ! End of the current iteration
      call Info_Mod_Iter_Print()

      if(ini >= min_ini) then
        call Control_Mod_Tolerance_For_Simple_Algorithm(simple_tol)
        if( u  % res <= simple_tol .and.  &
            v  % res <= simple_tol .and.  &
            w  % res <= simple_tol .and.  &
            mass_res <= simple_tol ) goto 1
      end if

    end do

    ! End of the current time step
1   call Info_Mod_Bulk_Print()

    ! Write the values in monitoring points
    if(.not. heat_transfer) then
      call Monitor_Mod_Write_4_Vars(n, u, v, w, p)
    else
      call Monitor_Mod_Write_5_Vars(n, u, v, w, t, p)
    end if

    ! Calculate mean values
    call Calculate_Mean(grid, n_stat, n)

    call User_Mod_Calculate_Mean(grid, n_stat, n)

    !-----------------------------------------------------!
    !   Recalculate the pressure drop                     !
    !   to keep the constant mass flux                    !
    !                                                     !
    !   First Newtons law:                                !
    !                                                     !
    !   F = m * a                                         !
    !                                                     !
    !   where:                                            !
    !                                                     !
    !   a = dv / dt = dFlux / dt * 1 / (A * rho)          !
    !   m = rho * v                                       !
    !   F = Pdrop * l * A = Pdrop * v                     !
    !                                                     !
    !   finally:                                          !
    !                                                     !
    !   Pdrop * v = rho * v * dFlux / dt * 1 / (A * rho)  !
    !                                                     !
    !   after cancelling: v and rho, it yields:           !
    !                                                     !
    !   Pdrop = dFlux/dt/A                                !
    !-----------------------------------------------------!
    if( abs(bulk % flux_x_o) >= TINY ) then
      bulk % p_drop_x = (bulk % flux_x_o - bulk % flux_x)  &
                         / (dt * bulk % area_x + TINY)
    end if
    if( abs(bulk % flux_y_o) >= TINY ) then
      bulk % p_drop_y = (bulk % flux_y_o - bulk % flux_y)  &
                         / (dt * bulk % area_y + TINY)
    end if
    if( abs(bulk % flux_z_o) >= TINY ) then
      bulk % p_drop_z = (bulk % flux_z_o - bulk % flux_z)  &
                         / (dt * bulk % area_z + TINY)
    end if

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
      call Backup_Mod_Save(grid, n, n_stat, name_save)
    end if

    ! Is it time to save results for post-processing
    if(save_now .or. exit_now .or. mod(n, rsi) .eq. 0) then
      call Comm_Mod_Wait
      call Save_Results(grid, name_save)

      ! Write results in user-customized format
      call User_Mod_Save_Results(grid, name_save)
    end if
   
    ! Just before the end of time step
    call User_Mod_End_Of_Time_Step(grid, n, time)

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
  call Save_Results(grid, name_save)
  call Backup_Mod_Save(grid, n, n_stat, name_save)

  ! Write results in user-customized format
  call User_Mod_Save_Results(grid, name_save)
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

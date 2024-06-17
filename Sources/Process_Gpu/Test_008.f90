#include "../Shared/Browse.h90"

!==============================================================================!
  subroutine Test_008
!------------------------------------------------------------------------------!
!>  Tests turbulent flow simulations
!------------------------------------------------------------------------------!
  use Native_Mod
  use Read_Controls_Mod
  use Process_Mod
  use Time_Mod
  use Iter_Mod
  use Turb_Mod
  use Gpu_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(Grid_Type)          :: Grid(MD)      ! computational grid
  type(Field_Type), target :: Flow(MD)      ! flow field
  type(Turb_Type)          :: Turb(MD)      ! turbulence models for flows
  real                     :: ts, te
  integer                  :: nc, ldt
  character(7)             :: root_control = 'control'
!==============================================================================!

  ! Start the parallel run and the profiler
  call Global % Start_Parallel
  call Profiler % Start('Test_008')

  O_Print '(a)', ' #====================================================='
  O_Print '(a)', ' # TEST 8: Turbulent flow simulations'
  O_Print '(a)', ' #====================================================='

  O_Print '(a)', ' # Opening the control file '//root_control
  call Control % Open_Root_File(root_control)
  call Control % Time_Step(Flow(1) % dt, verbose=.true.)

  ! Initialize Info_Mod
  call Info % Start_Info()

  O_Print '(a)', ' # Creating a grid'
  call Grid(1) % Load_And_Prepare_For_Processing(1)

  O_Print '(a)', ' # Reading physical models'
  call Read_Control % Physical_Models(Flow(1), Turb(1))

  nc = Grid(1) % n_cells
  O_Print '(a, i12)',   ' # The problem size is: ', nc
  O_Print '(a,es12.3)', ' # Solver tolerace is : ', PICO

  O_Print '(a)', ' #----------------------------------------------------'
  O_Print '(a)', ' # Be careful with memory usage.  If you exceed the'
  O_Print '(a)', ' # 90% (as a rule of thumb) of the memory your GPU'
  O_Print '(a)', ' # card has the program will become memory bound no'
  O_Print '(a)', ' # matter how you wrote it, and it may even crash.'
  O_Print '(a)', ' #----------------------------------------------------'

  O_Print '(a)', ' # Creating a flow field'
  call Flow(1) % Create_Field(Grid(1))

  O_Print '(a)', ' # Creating turbulent models space'
  call Turb(1)  % Create_Turb(Flow(1), Grid(1))

  O_Print '(a)', ' # Reading physical properties'
  call Read_Control % Physical_Properties(Flow(1), Grid(1))

  ! I am not sure when to call this, but this is a good guess
  call Read_Control % Boundary_Conditions(Flow(1), Grid(1))

  ! Read numerical models from control file (after the memory is allocated)
  call Read_Control % Numerical_Schemes(Flow(1), Grid(1))

  O_Print '(a)', ' # Reading native solvers'
  call Read_Control % Native_Solvers(Flow(1), Grid(1))

  O_Print '(a)', ' # Calculating gradient matrix for the field'
  call Flow(1) % Calculate_Grad_Matrix(Grid(1))

  ! Initialize variables
  call Process % Initialize_Variables(Turb(1), Flow(1), Grid(1))

  !----------------------------------------------------------!
  !   Copy all useful data to the device, that means grid,   !
  !   field and solvers                                      !
  !----------------------------------------------------------!
  call Gpu % Grid_Copy_To_Device(Turb(1), Grid(1))
  call Gpu % Field_Copy_To_Device(Flow(1))
  call Gpu % Native_Copy_To_Device(Flow(1) % Nat)
  call Gpu % Turb_Copy_To_Device(Turb(1), Flow(1))

  !------------------------------------------!
  !                                          !
  !   Performing a time loop on the device   !
  !                                          !
  !------------------------------------------!

  ! Read the number of time steps ...
  call Control % Number_Of_Time_Steps(ldt, verbose=.true.)

  ! ... and if iterations from control file
  call Read_Control % Iterations()

  ! Time stepping initializations
  call Time % Set_Curr_Dt(0)
  call Time % Set_First_Dt(0)
  call Time % Set_Last_Dt(ldt)

  call Control % Results_Save_Interval     (Results % interval, verbose=.true.)
  call Control % Save_Results_At_Boundaries(Results % boundary)

  ! Allocate CPU memory for working arrays (currently used for saving)
  call Work % Allocate_Work(Grid, n_r_cell=12,  n_r_face=0,  n_r_node=0,  &
                                  n_i_cell= 6,  n_i_face=0,  n_i_node=0)

  O_Print '(a)', ' # Performing a demo of the computing momentum equations'
  call cpu_time(ts)
  do while (Time % Needs_More_Steps())

    ! Start info boxes
    call Info % Time_Start()
    call Info % Iter_Start()
    call Info % Bulk_Start()

    call Info % Time_Fill(Time % Curr_Dt(), Time % Get_Time())
    call Info % Time_Print()

    ! Turbulence models initializations
    call Turb(1) % Init_Turb(Flow(1), Grid(1))

    !-----------------------------------!
    !   Iterations within a time step   !
    !-----------------------------------!
    do while (Iter % Needs_More_Iterations(Flow, 1))

      ! Beginning of iteration
      call Info % Iter_Fill(Iter % Current())

      ! New in GPU version: compute pressure gradients here
      call Flow(1) % Grad_Pressure(Grid(1), Flow(1) % p)

      call Process % Compute_Momentum(Flow(1), Grid(1), comp=1)
      call Process % Compute_Momentum(Flow(1), Grid(1), comp=2)
      call Process % Compute_Momentum(Flow(1), Grid(1), comp=3)

      call Process % Compute_Pressure(Flow(1), Grid(1))

      ! Correct velocity components
      call Flow(1) % Grad_Pressure(Grid(1), Flow(1) % pp)
      call Process % Correct_Velocity(Flow(1), Grid(1))

      ! Deal with turbulence
      call Turb(1) % Main_Turb(Flow(1), Grid(1))

      if(Flow(1) % heat_transfer) then
        call Process % Compute_Energy(Flow(1), Grid(1))
      end if

      call Process % Update_Boundary_Values(Flow(1), Grid(1), 'ALL')

      ! End of the current iteration
      call Info % Iter_Print(1)

    end do  ! iterations

    ! Calculate bulk fluxes and adjust pressure drops
    call Flow(1) % Calculate_Bulk_Velocities(Grid(1))
    call Flow(1) % Adjust_P_Drops(Grid(1))

    ! Print the bulk values from the Info_Mod
    call Info % Bulk_Print(Flow(1), 1, 1)

    if(mod(Time % Curr_Dt(), Results % interval) .eq. 0) then
      call Gpu % Turb_Update_Host(Turb(1), Flow(1))
      call Gpu % Field_Update_Host(Flow(1))
      call Gpu % Grid_Update_Host(Turb(1), Grid(1))
      call Results % Main_Results(Turb(1), Flow(1), Grid(1), 1)
    end if

  end do    ! time steps
  call cpu_time(te)

  ! Save results
  call Gpu % Turb_Update_Host(Turb(1), Flow(1))
  call Gpu % Field_Update_Host(Flow(1))
  call Gpu % Grid_Update_Host(Turb(1), Grid(1))
  call Results % Main_Results(Turb(1), Flow(1), Grid(1), 1)

  call Work % Finalize_Work()

  ! End the profiler and the parallel run
  call Profiler % Stop('Test_008')
  call Profiler % Statistics(indent=24)
  call Global % End_Parallel

  end subroutine

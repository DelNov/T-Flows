#include "../Shared/Browse.h90"

!==============================================================================!
  subroutine Test_009
!------------------------------------------------------------------------------!
!>  Tests calling of the CG algorithm from the Native_Mod
!------------------------------------------------------------------------------!
  use Native_Mod
  use Read_Controls_Mod
  use Process_Mod
  use Time_Mod
  use Iter_Mod
  use Gpu_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(Grid_Type)          :: Grid(MD)      ! computational grid
  type(Field_Type), target :: Flow(MD)      ! flow field
  real                     :: ts, te
  integer                  :: nc, c, ldt
  logical                  :: exit_now
  character(11)            :: name_vel     = 'TTTT_II_uvw'
  character( 9)            :: name_p       = 'TTTT_II_p'
!@character(10)            :: name_pp      = 'TTTT_II_pp'
!@character(14)            :: name_grad_p  = 'TTTT_II_grad_p'
!@character(15)            :: name_grad_pp = 'TTTT_II_grad_pp'
  character(11)            :: root_control = 'control.009'
!==============================================================================!

  ! Start the parallel run and the profiler
  call Global % Start_Parallel
  call Profiler % Start('Test_009')

  O_Print '(a)', ' #====================================================='
  O_Print '(a)', ' # TEST 9: Flow over a backward facing step'
  O_Print '(a)', ' #====================================================='

  O_Print '(a)', ' # Opening the control file '//root_control
  call Control % Open_Root_File(root_control)
  call Control % Time_Step(Flow(1) % dt, verbose=.true.)

  ! Initialize Info_Mod
  call Info % Start_Info()

  O_Print '(a)', ' # Creating a grid'
  call Grid(1) % Load_And_Prepare_For_Processing(1)

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

  O_Print '(a)', ' # Reading physical properties'
  call Read_Control % Physical_Properties(Flow(1))

  ! I am not sure when to call this, but this is a good guess
  call Read_Control % Boundary_Conditions(Flow(1))

  ! Read numerical models from control file (after the memory is allocated)
  call Read_Control % Numerical_Schemes(Flow(1))

  O_Print '(a)', ' # Reading native solvers'
  call Read_Control % Native_Solvers(Flow(1))

  ! Form preconditioning matrices on host
  ! (Must be before transferring them)
  call Flow(1) % Nat % Prec_Form(Flow(1) % Nat % C, Flow(1) % Nat % M)
  call Flow(1) % Nat % Prec_Form(Flow(1) % Nat % C, Flow(1) % Nat % A)

  O_Print '(a)', ' # Calculating gradient matrix for the field'
  call Flow(1) % Calculate_Grad_Matrix()

  ! You are going to need connectivity matrix on device
  call Gpu % Sparse_Con_Copy_To_Device(Flow(1) % Nat % C)

  ! OK, once you formed the preconditioners, you
  ! will want to keep these matrices on the device
  call Gpu % Sparse_Val_Copy_To_Device(Flow(1) % Nat % M)
  call Gpu % Sparse_Val_Copy_To_Device(Flow(1) % Nat % A)

  ! and that bloody right-hand-side vector too
  call Gpu % Vector_Real_Copy_To_Device(Flow(1) % Nat % b)

  ! In addition to system matrices of your discretized
  ! equations, you will want to have gradient matrices, as
  ! well as cell connectivity and cell coordinates on the
  ! device (they are all needed for gradients), ...
  call Gpu % Matrix_Real_Copy_To_Device(Flow(1) % grad_c2c)
  call Gpu % Matrix_Int_Copy_To_Device(Grid(1) % faces_c)
  call Gpu % Vector_Int_Copy_To_Device(Grid(1) % cells_n_cells)
  call Gpu % Matrix_Int_Copy_To_Device(Grid(1) % cells_c)
  call Gpu % Matrix_Int_Copy_To_Device(Grid(1) % cells_f)
  call Gpu % Vector_Real_Copy_To_Device(Grid(1) % xc)
  call Gpu % Vector_Real_Copy_To_Device(Grid(1) % yc)
  call Gpu % Vector_Real_Copy_To_Device(Grid(1) % zc)
  call Gpu % Vector_Real_Copy_To_Device(Grid(1) % dx)
  call Gpu % Vector_Real_Copy_To_Device(Grid(1) % dy)
  call Gpu % Vector_Real_Copy_To_Device(Grid(1) % dz)
  call Gpu % Vector_Real_Copy_To_Device(Grid(1) % sx)
  call Gpu % Vector_Real_Copy_To_Device(Grid(1) % sy)
  call Gpu % Vector_Real_Copy_To_Device(Grid(1) % sz)
  call Gpu % Vector_Real_Copy_To_Device(Grid(1) % s)
  call Gpu % Vector_Real_Copy_To_Device(Grid(1) % vol)
  call Gpu % Vector_Int_Copy_To_Device(Grid(1) % region % f_face)
  call Gpu % Vector_Int_Copy_To_Device(Grid(1) % region % l_face)
  call Gpu % Vector_Int_Copy_To_Device(Grid(1) % region % f_cell)
  call Gpu % Vector_Int_Copy_To_Device(Grid(1) % region % l_cell)

  ! ... and the vectors of the native suite of solvers
  call Gpu % Native_Transfer_To_Device(Flow(1) % Nat)

  ! OK, fine, now you have all sort of matrices and supporting
  ! data on the device, but you will also need variables sol-
  ! ved for (pp % n, u % n, v % n and w % n), universal source
  ! for them all (b) and the variables whose gradients are
  ! being computed (pp % n and p % n) as well as gradient com
  ! ponents (pp % x, pp % y, pp % z, p % x, p % y and p % z)
  call Gpu % Vector_Real_Copy_To_Device(Flow(1) % pp % n)
  call Gpu % Vector_Real_Copy_To_Device(Flow(1) % p % n)
  call Gpu % Vector_Real_Copy_To_Device(Flow(1) % u % n)
  call Gpu % Vector_Real_Copy_To_Device(Flow(1) % v % n)
  call Gpu % Vector_Real_Copy_To_Device(Flow(1) % w % n)
  call Gpu % Vector_Real_Copy_To_Device(Flow(1) % u % o)
  call Gpu % Vector_Real_Copy_To_Device(Flow(1) % v % o)
  call Gpu % Vector_Real_Copy_To_Device(Flow(1) % w % o)
  call Gpu % Vector_Real_Copy_To_Device(Flow(1) % pp % x)
  call Gpu % Vector_Real_Copy_To_Device(Flow(1) % pp % y)
  call Gpu % Vector_Real_Copy_To_Device(Flow(1) % pp % z)
  call Gpu % Vector_Real_Copy_To_Device(Flow(1) % p % x)
  call Gpu % Vector_Real_Copy_To_Device(Flow(1) % p % y)
  call Gpu % Vector_Real_Copy_To_Device(Flow(1) % p % z)
  call Gpu % Vector_Real_Copy_To_Device(Flow(1) % v_flux)

  ! This should be done for each domain, whenever a new domain is solved
  call Flow(1) % Update_Aliases()

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
  call Work % Allocate_Work(Grid, n_r_cell=8,  n_r_face=0,  n_r_node=0,  &
                                  n_i_cell=6,  n_i_face=0,  n_i_node=0)

  O_Print '(a)', ' # Performing a demo of the computing momentum equations'
  call cpu_time(ts)
  do while (Time % Needs_More_Steps())

    ! Start info boxes
    call Info % Time_Start()
    call Info % Iter_Start()
    call Info % Bulk_Start()

    call Info % Time_Fill(Time % Curr_Dt(), Time % Get_Time())
    call Info % Time_Print()

    ! Preparation for the new time step
    !$acc parallel loop
    do c = 1, nc
      u_o(c) = u_n(c)
      v_o(c) = v_n(c)
      w_o(c) = w_n(c)
    end do
    !$acc end parallel

    write(name_vel    (1:4), '(i4.4)') Time % Curr_Dt()
    write(name_p      (1:4), '(i4.4)') Time % Curr_Dt()
!@  write(name_pp     (1:4), '(i4.4)') Time % Curr_Dt()
!@  write(name_grad_pp(1:4), '(i4.4)') Time % Curr_Dt()
!@  write(name_grad_p (1:4), '(i4.4)') Time % Curr_Dt()

    !-----------------------------------!
    !   Iterations within a time step   !
    !-----------------------------------!
    do while (Iter % Needs_More_Iterations(Flow, 1))

      ! Beginning of iteration
      call Info % Iter_Fill(Iter % Current())

!@    write(name_vel    (6:7), '(i2.2)') Iter % Current()
!@    write(name_p      (6:7), '(i2.2)') Iter % Current()
!@    write(name_pp     (6:7), '(i2.2)') Iter % Current()
!@    write(name_grad_pp(6:7), '(i2.2)') Iter % Current()
!@    write(name_grad_p (6:7), '(i2.2)') Iter % Current()

      call Process % Compute_Momentum(Flow(1), Grid(1), comp=1)
      call Process % Compute_Momentum(Flow(1), Grid(1), comp=2)
      call Process % Compute_Momentum(Flow(1), Grid(1), comp=3)

      call Process % Compute_Pressure(Flow(1), Grid(1))

      call Flow(1) % Grad_Pressure(Grid(1), Flow(1) % pp)
      call Process % Correct_Velocity(Flow(1), Grid(1))
      call Flow(1) % Grad_Pressure(Grid(1), Flow(1) % p)

      call Process % Update_Boundary_Values(Flow(1), Grid(1), 'ALL')

      ! End of the current iteration
      call Info % Iter_Print(1)

    end do  ! iterations


    ! Print the bulk values from the Info_Mod
    call Info % Bulk_Print(Flow(1), 1, 1)

    if(mod(Time % Curr_Dt(), Results % interval) .eq. 0) then
      call Gpu % Vector_Update_Host(Flow(1) % u % n)
      call Gpu % Vector_Update_Host(Flow(1) % v % n)
      call Gpu % Vector_Update_Host(Flow(1) % w % n)
      call Gpu % Vector_Update_Host(Flow(1) % p % n)
      call Results % Main_Results(1, Flow(1), exit_now)
    end if

  end do    ! time steps
  call cpu_time(te)

  ! Save results
  call Gpu % Vector_Update_Host(Flow(1) % u % n)
  call Gpu % Vector_Update_Host(Flow(1) % v % n)
  call Gpu % Vector_Update_Host(Flow(1) % w % n)
  call Gpu % Vector_Update_Host(Flow(1) % p % n)
  call Results % Main_Results(1, Flow(1), exit_now)

  call Work % Finalize_Work()

  ! End the profiler and the parallel run
  call Profiler % Stop('Test_009')
  call Profiler % Statistics(indent=24)
  call Global % End_Parallel

  end subroutine

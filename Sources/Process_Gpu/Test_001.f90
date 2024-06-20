#include "../Shared/Browse.h90"
#include "../Shared/Macros.h90"

!==============================================================================!
  subroutine Test_001
!------------------------------------------------------------------------------!
!>  Tests matrix-vector product
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Linalg_Mod
  use Process_Mod
  use Gpu_Mod
  use Read_Controls_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(Sparse_Con_Type), pointer :: Acon
  type(Sparse_Val_Type), pointer :: Aval
  real, allocatable              :: b(:), c(:)
  real,      contiguous, pointer :: temp(:)
  type(Grid_Type)                :: Grid
  type(Field_Type),       target :: Flow            ! flow field
  integer,             parameter :: N_STEPS = 1200  ! spend some time on device
  integer                        :: nc, ni, step
  character(len=11)              :: root_control = 'control.001'
!==============================================================================!

  ! Start the parallel run and the profiler
  call Global % Start_Parallel
  call Profiler % Start('Test_001')

  call Work % Connect_Real_Cell(temp)

  O_Print '(a)', ' #===================================================='
  O_Print '(a)', ' # TEST 1: Performing a sparse-matrix vector product'
  O_Print '(a)', ' #===================================================='

  O_Print '(a)', ' # Opening the control file '//root_control
  call Control % Open_Root_File(root_control)

  O_Print '(a)', ' # Creating a grid'
  call Grid % Load_And_Prepare_For_Processing(1)

  nc = Grid % n_cells
  ni = Grid % n_cells - Grid % Comm % n_buff_cells
  O_Print '(a,i12)', ' # The problem size is: ', ni

  O_Print '(a)', ' #----------------------------------------------------'
  O_Print '(a)', ' # Be careful with memory usage.  If you exceed the'
  O_Print '(a)', ' # 90% (as a rule of thumb) of the memory your GPU'
  O_Print '(a)', ' # card has the program will become memory bound no'
  O_Print '(a)', ' # matter how you wrote it, and it may even crash.'
  O_Print '(a)', ' #----------------------------------------------------'

  O_Print '(a)', ' # Creating a grid'
  call Flow % Create_Field(Grid)

  O_Print '(a)', ' # Reading physical properties'
  call Read_Control % Physical_Properties(Flow, Grid)

  ! Discretize the matrix for diffusion
  temp(:) = 1.0
  call Process % Form_System_Matrix(Acon, Aval, Flow, Grid,  &
                                    Flow % density, temp,    &
                                    Flow % viscosity,        &
                                    1.0)

  ! Take the alias now
  Acon => Flow % Nat % C
  Aval => Flow % Nat % A(MATRIX_UVW)

  allocate(b(nc))
  allocate(c(nc))

  ! Fill up the right-hand side vector up to buffers
  b(1:nc) = 0.0
  b(1:ni) = 2.0

  ! Copy operand matrix and vector to the device ...
  ! ... and reserve memory for result vector on device
  call Gpu % Sparse_Con_Copy_To_Device(Acon)
  call Gpu % Sparse_Val_Copy_To_Device(Aval)
  call Gpu % Vector_Real_Copy_To_Device(b)
  call Gpu % Vector_Real_Create_On_Device(c)

  !-----------------------------------------------!
  !   Performing a fake time loop on the device   !
  !-----------------------------------------------!
  O_Print '(a,i6,a)', ' # Performing ', N_STEPS, ' sparse-matrix vector products'
  call Profiler % Start('Useful_Work')
  do step = 1, N_STEPS
    call Linalg % Mat_X_Vec(ni, c(1:nc), Acon, Aval, b(1:nc))
  end do
  call Profiler % Stop('Useful_Work')

  ! Mat_X_Vec gives good values in inside cells, but not in buffers. For
  ! linear solvers those might be sufficient for the rest of the algorithm,
  ! but here we want to plot results, so refreshing the buffer is needed
  call Grid % Exchange_Inside_Cells_Real(c(1:nc))

  ! Copy results back to host
  call Gpu % Vector_Update_Host(c)

  ! Destroy data on the device, you don't need them anymore
  call Gpu % Sparse_Con_Destroy_On_Device(Acon)
  call Gpu % Sparse_Val_Destroy_On_Device(Aval)
  call Gpu % Vector_Real_Destroy_On_Device(b)
  call Gpu % Vector_Real_Destroy_On_Device(c)

  ! Print result
  O_Print '(a,es12.3)', ' vector c(1  ):', c(1)
  O_Print '(a,es12.3)', ' vector c(2  ):', c(2)
  O_Print '(a,es12.3)', ' vector c(n-1):', c(Grid % n_cells - 1)
  O_Print '(a,es12.3)', ' vector c(n  ):', c(Grid % n_cells)

  ! Save results
  call Grid % Save_Debug_Vtu("result",              &
                             inside_name="Result",  &
                             inside_cell=c)

  call Work % Disconnect_Real_Cell(temp)

  ! End the profiler and the parallel run
  call Profiler % Stop('Test_001')
  call Profiler % Statistics(indent=24)
  call Global % End_Parallel

  end subroutine

#include "../Shared/Browse.h90"

!==============================================================================!
  subroutine Test_004
!------------------------------------------------------------------------------!
!>  Tests calling of the CG algorithm from the Native_Mod
!------------------------------------------------------------------------------!
  use Native_Mod
  use Read_Controls_Mod
  use Process_Mod
  use Gpu_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(Sparse_Con_Type), pointer :: Acon
  type(Sparse_Val_Type), pointer :: Aval
  type(Grid_Type)                :: Grid  ! computational grid
  type(Field_Type),       target :: Flow  ! flow field
  real,                  pointer :: b(:), x(:)
  integer                        :: nc, n
  real                           :: fin_res
  character(len=11)              :: root_control = 'control.004'
!==============================================================================!

  ! Start the parallel run and the profiler
  call Global % Start_Parallel
  call Profiler % Start('Test_004')

  O_Print '(a)', ' #====================================================='
  O_Print '(a)', ' # TEST 4: Call Conjugate Gradient from Native_Mod'
  O_Print '(a)', ' #====================================================='

  O_Print '(a)', ' # Opening the control file '//root_control
  call Control % Open_Root_File(root_control)

  O_Print '(a)', ' # Creating a grid'
  call Grid % Load_And_Prepare_For_Processing(1)

  nc = Grid % n_cells
  O_Print '(a,i12)',    ' # The problem size is: ', nc
  O_Print '(a,es12.3)', ' # Solver tolerace is : ', PICO

  O_Print '(a)', ' #----------------------------------------------------'
  O_Print '(a)', ' # Be careful with memory usage.  If you exceed the'
  O_Print '(a)', ' # 90% (as a rule of thumb) of the memory your GPU'
  O_Print '(a)', ' # card has the program will become memory bound no'
  O_Print '(a)', ' # matter how you wrote it, and it may even crash.'
  O_Print '(a)', ' #----------------------------------------------------'

  O_Print '(a)', ' # Creating a field'
  call Flow % Create_Field(Grid)

  O_Print '(a)', ' # Reading physical properties'
  call Read_Control % Physical_Properties(Flow)

  ! I am not sure when to call this, but this is a good guess
  call Read_Control % Boundary_Conditions(Flow)

  ! Read numerical models from control file (after the memory is allocated)
  call Read_Control % Numerical_Schemes(Flow)

  ! Transfer the necessary grid components to the device
  call Gpu % Matrix_Int_Copy_To_Device(Grid % faces_c)
  call Gpu % Vector_Int_Copy_To_Device(Grid % cells_n_cells)
  call Gpu % Matrix_Int_Copy_To_Device(Grid % cells_c)
  call Gpu % Matrix_Int_Copy_To_Device(Grid % cells_f)
  call Gpu % Vector_Int_Copy_To_Device(Grid % region % f_face)
  call Gpu % Vector_Int_Copy_To_Device(Grid % region % l_face)
  call Gpu % Vector_Int_Copy_To_Device(Grid % region % f_cell)
  call Gpu % Vector_Int_Copy_To_Device(Grid % region % l_cell)

  ! Important before calling functions in Process which are ported to GPUs)
  call Flow % Update_Aliases()

  ! Take the aliases now
  Acon => Flow % Nat % C
  Aval => Flow % Nat % A(MATRIX_UVW)
  b    => Flow % Nat % b
  x    => Flow % u % n

  ! Initialize solution (start from 1 not to overwrite boundary conditions)
  x(1:nc) = 0.0

  ! Copy components of the linear system to the device
  call Gpu % Sparse_Con_Copy_To_Device(Acon)
  call Gpu % Sparse_Val_Copy_To_Device(Aval)
  call Gpu % Vector_Real_Copy_To_Device(x)
  call Gpu % Vector_Real_Copy_To_Device(b)

  ! Allocate vectors related to CG algorithm on the device
  call Gpu % Native_Transfer_To_Device(Flow % Nat)

  ! Copy physical properties as well
  call Gpu % Vector_Real_Copy_To_Device(Flow % viscosity)
  call Gpu % Vector_Real_Copy_To_Device(Flow % density)

  !-------------------------------------------------!
  !   Discretize the linear system for conduction   !
  !-------------------------------------------------!
  call Process % Form_System_Matrix(Acon, Aval, Flow, Grid,       &
                                    Flow % density, Flow % ones,  &
                                    Flow % viscosity)
  call Process % Insert_Momentum_Bc(Flow, Grid, comp=1)

  !-----------------------------------------------!
  !   Performing a fake time loop on the device   !
  !-----------------------------------------------!
  O_Print '(a)', ' # Performing a demo of the preconditioned CG method'
  call Profiler % Start('Useful_Work')
  call Flow % Nat % Cg(Acon, Aval, x, b, nc, n, PICO, fin_res)
  call Profiler % Stop('Useful_Work')

  ! Copy results back to host
  call Gpu % Vector_Update_Host(x)

  ! Destroy all data on the device, you don't need them anymore
  call Gpu % Matrix_Int_Destroy_On_Device(Grid % faces_c)
  call Gpu % Vector_Int_Destroy_On_Device(Grid % cells_n_cells)
  call Gpu % Matrix_Int_Destroy_On_Device(Grid % cells_c)
  call Gpu % Matrix_Int_Destroy_On_Device(Grid % cells_f)
  call Gpu % Vector_Int_Destroy_On_Device(Grid % region % f_face)
  call Gpu % Vector_Int_Destroy_On_Device(Grid % region % l_face)
  call Gpu % Vector_Int_Destroy_On_Device(Grid % region % f_cell)
  call Gpu % Vector_Int_Destroy_On_Device(Grid % region % l_cell)

  call Gpu % Sparse_Con_Destroy_On_Device(Acon)
  call Gpu % Sparse_Val_Destroy_On_Device(Aval)
  call Gpu % Vector_Real_Destroy_On_Device(x)
  call Gpu % Vector_Real_Destroy_On_Device(b)

  call Gpu % Native_Destroy_On_Device(Flow % Nat)

  ! Print result
  O_Print '(a,es12.3)', ' vector u(1  ):', x(1)
  O_Print '(a,es12.3)', ' vector u(2  ):', x(2)
  O_Print '(a,es12.3)', ' vector u(3  ):', x(3)
  O_Print '(a,es12.3)', ' vector u(n-2):', x(Grid % n_cells-2)
  O_Print '(a,es12.3)', ' vector u(n-1):', x(Grid % n_cells-1)
  O_Print '(a,es12.3)', ' vector u(n  ):', x(Grid % n_cells)

  ! Save results
  call Grid % Save_Debug_Vtu("result",              &
                             scalar_name="Result",  &
                             scalar_cell=x)

  ! End the profiler and the parallel run
  call Profiler % Stop('Test_004')
  call Profiler % Statistics(indent=24)
  call Global % End_Parallel

  end subroutine

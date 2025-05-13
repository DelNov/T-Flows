#include "../../Shared/Browse.h90"
#include "../../Shared/Macros.h90"

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
  type(Grid_Type)                :: Grid(1)  ! computational grid
  type(Field_Type),       target :: Flow     ! flow field
  type(Turb_Type), target        :: Turb
  real,      contiguous, pointer :: b(:), x(:)
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
  call Grid(1) % Load_And_Prepare_For_Processing(1)
  call Grid(1) % Copy_Grid_To_Device()

  nc = Grid(1) % n_cells
  O_Print '(a,i12)',    ' # The problem size is: ', nc
  O_Print '(a,es12.3)', ' # Solver tolerace is : ', PICO

# if T_FLOWS_GPU == 1
    O_Print '(a)', ' #----------------------------------------------------'
    O_Print '(a)', ' # Be careful with memory usage.  If you exceed the'
    O_Print '(a)', ' # 90% (as a rule of thumb) of the memory your GPU'
    O_Print '(a)', ' # card has the program will become memory bound no'
    O_Print '(a)', ' # matter how you wrote it, and it may even crash.'
    O_Print '(a)', ' #----------------------------------------------------'
# endif

  O_Print '(a)', ' # Creating a field'
  call Flow % Create_Field(Grid(1))

  O_Print '(a)', ' # Reading physical properties'
  call Read_Control % Physical_Properties(Grid(1), Flow)

  ! I am not sure when to call this, but this is a good guess
  call Read_Control % Boundary_Conditions(Grid(1), Flow, Turb)

  ! Read numerical models from control file (after the memory is allocated)
  call Read_Control % Numerical_Schemes(Grid(1), Flow, Turb)

  ! Take the aliases now
  Acon => Flow % Nat % C
  Aval => Flow % Nat % A
  b    => Flow % Nat % b
  x    => Flow % u % n

  ! Initialize solution (start from 1 not to overwrite boundary conditions)
  x(1:nc) = 0.0

  ! Copy components of the linear system to the device
  call Acon % Copy_Sparse_Con_To_Device()
  call Aval % Copy_Sparse_Val_To_Device()
  call Gpu % Vector_Real_Copy_To_Device(x)
  call Gpu % Vector_Real_Copy_To_Device(b)

  ! Allocate vectors related to CG algorithm on the device
  call Flow % Nat % Copy_Native_To_Device()

  ! Allocate CPU memory for working arrays (currently used for saving)
  call Work % Allocate_Work(Grid, n_r_cell=24,  n_r_face=0,  n_r_node=0,  &
                                  n_i_cell= 6,  n_i_face=0,  n_i_node=0)
  call Work % Create_Work_On_Device()

  ! Copy physical properties as well
  call Flow % Copy_Field_To_Device()

  !-------------------------------------------------!
  !   Discretize the linear system for conduction   !
  !-------------------------------------------------!
  call Process % Form_Momentum_Matrix(Grid(1), Flow, Turb, Flow % viscosity, 1.0)
  call Process % Insert_Momentum_Bc(Grid(1), Flow, comp=1)

  !-----------------------------------------------!
  !   Performing a fake time loop on the device   !
  !-----------------------------------------------!
  O_Print '(a)', ' # Performing a demo of the preconditioned CG method'
  call Profiler % Start('Useful_Work')
  call Flow % Nat % Cg(x, nc, n, PICO, fin_res)
  call Profiler % Stop('Useful_Work')

  ! Copy results back to host
  call Gpu % Vector_Update_Host(x)

  ! Destroy all data on the device, you don't need them anymore
  call Acon % Destroy_Sparse_Con_On_Device()
  call Aval % Destroy_Sparse_Val_On_Device()
  call Gpu % Vector_Real_Destroy_On_Device(x)
  call Gpu % Vector_Real_Destroy_On_Device(b)

  call Flow % Nat % Destroy_Native_On_Device()

  ! Print result
  O_Print '(a,es12.3)', ' vector u(1  ):', x(1)
  O_Print '(a,es12.3)', ' vector u(2  ):', x(2)
  O_Print '(a,es12.3)', ' vector u(3  ):', x(3)
  O_Print '(a,es12.3)', ' vector u(n-2):', x(Grid(1) % n_cells-2)
  O_Print '(a,es12.3)', ' vector u(n-1):', x(Grid(1) % n_cells-1)
  O_Print '(a,es12.3)', ' vector u(n  ):', x(Grid(1) % n_cells)

  ! Save results
  call Grid(1) % Save_Debug_Vtu("result",              &
                                scalar_name="Result",  &
                                scalar_cell=x)

  ! End the profiler and the parallel run
  call Profiler % Stop('Test_004')
  call Profiler % Statistics(indent=24)
  call Global % End_Parallel

  end subroutine

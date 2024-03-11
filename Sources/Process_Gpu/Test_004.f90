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
  real                           :: ts, te
  integer                        :: n
  character(len=11)              :: root_control = 'control.004'
!==============================================================================!

  print '(a)', ' #====================================================='
  print '(a)', ' # TEST 4: Call Conjugate Gradient from Native_Mod'
  print '(a)', ' #====================================================='

  print '(a)', ' # Opening the control file '//root_control
  call Control % Open_Root_File(root_control)

  print '(a)', ' # Creating a grid'
  call Grid % Load_And_Prepare_For_Processing(1)

  n = Grid % n_cells
  print '(a,i12)',    ' # The problem size is: ', n
  print '(a,es12.3)', ' # Solver tolerace is : ', PICO

  print '(a)', ' #----------------------------------------------------'
  print '(a)', ' # Be careful with memory usage.  If you exceed the'
  print '(a)', ' # 90% (as a rule of thumb) of the memory your GPU'
  print '(a)', ' # card has the program will become memory bound no'
  print '(a)', ' # matter how you wrote it, and it may even crash.'
  print '(a)', ' #----------------------------------------------------'

  print '(a)', ' # Creating a field'
  call Flow % Create_Field(Grid)

  ! I am not sure when to call this, but this is a good guess
  call Read_Control % Boundary_Conditions(Flow)

  ! Discretize the matrix for diffusion
  call Process % Form_Diffusion_Matrix(Flow)

  call Gpu % Matrix_Int_Copy_To_Device(Grid % faces_c)
  call Gpu % Vector_Real_Copy_To_Device(Grid % s)
  call Gpu % Vector_Real_Copy_To_Device(Grid % d)
  call Gpu % Vector_Int_Copy_To_Device(Grid % region % f_face)
  call Gpu % Vector_Int_Copy_To_Device(Grid % region % l_face)

  ! Important before calling functions in Process which are ported to GPUs)
  call Flow % Update_Aliases()

  call Process % Insert_Diffusion_Bc(Flow, comp=1)

  ! Take the aliases now
  Acon => Flow % Nat % C
  Aval => Flow % Nat % M
  b    => Flow % Nat % b
  x    => Flow % u % n

  ! Initialize solution
  x(:) = 0.0

  ! Before copying matrix components, create a preconditioning diagonal
  call Flow % Nat % Prec_Form(Acon, Aval)

  ! Copy components of the linear system to the device
  call Gpu % Sparse_Con_Copy_To_Device(Acon)
  call Gpu % Sparse_Val_Copy_To_Device(Aval)
  call Gpu % Vector_Real_Copy_To_Device(x)
  call Gpu % Vector_Real_Copy_To_Device(b)

  ! Allocate vectors related to CG algorithm on the device
  call Gpu % Native_Transfer_To_Device(Flow % Nat)

  !-----------------------------------------------!
  !   Performing a fake time loop on the device   !
  !-----------------------------------------------!
  print '(a)', ' # Performing a demo of the preconditioned CG method'
  call cpu_time(ts)
  call Flow % Nat % Cg(Acon, Aval, x, b, n, PICO)
  call cpu_time(te)

  ! Copy results back to host
  call Gpu % Vector_Update_Host(x)

  ! Destroy all data on the device, you don't need them anymore
  call Gpu % Matrix_Int_Destroy_On_Device(Grid % faces_c)
  call Gpu % Vector_Real_Destroy_On_Device(Grid % s)
  call Gpu % Vector_Real_Destroy_On_Device(Grid % d)
  call Gpu % Vector_Int_Destroy_On_Device(Grid % region % f_face)
  call Gpu % Vector_Int_Destroy_On_Device(Grid % region % l_face)

  call Gpu % Sparse_Con_Destroy_On_Device(Acon)
  call Gpu % Sparse_Val_Destroy_On_Device(Aval)
  call Gpu % Vector_Real_Destroy_On_Device(x)
  call Gpu % Vector_Real_Destroy_On_Device(b)

  call Gpu % Native_Destroy_On_Device(Flow % Nat)

  ! Print result
  print '(a,es12.3)', ' vector u(1  ):', x(1)
  print '(a,es12.3)', ' vector u(2  ):', x(2)
  print '(a,es12.3)', ' vector u(3  ):', x(3)
  print '(a,es12.3)', ' vector u(n-2):', x(Grid % n_cells-2)
  print '(a,es12.3)', ' vector u(n-1):', x(Grid % n_cells-1)
  print '(a,es12.3)', ' vector u(n  ):', x(Grid % n_cells)

  ! Save results
  call Grid % Save_Debug_Vtu("solution",       &
                             scalar_name="u",  &
                             scalar_cell=x)

  print '(a,f12.3,a)', ' # Time elapsed for TEST 4: ', te-ts, ' [s]'

  end subroutine

!==============================================================================!
  subroutine Test_001
!------------------------------------------------------------------------------!
!>  Tests matrix-vector product
!------------------------------------------------------------------------------!
  use Linalg_Mod
  use Process_Mod
  use Gpu_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(Sparse_Type), pointer :: A
  real, allocatable          :: b(:), c(:)
  type(Grid_Type)            :: Grid
  type(Field_Type),   target :: Flow                   ! flow field
  integer,         parameter :: N_STEPS = 1200  ! spend enough time on device
  real                       :: ts, te
  integer                    :: n, step
  character(len=11)          :: root_control    = 'control.001'
!==============================================================================!

  print '(a)', ' #===================================================='
  print '(a)', ' # TEST 1: Performing a sparse-matrix vector product'
  print '(a)', ' #===================================================='

  print '(a)', ' # Opening the control file '//root_control
  call Control % Open_Root_File(root_control)

  print '(a)', ' # Creating a grid'
  call Grid % Load_And_Prepare_For_Processing(1)

  n = Grid % n_cells
  print '(a,i12)', ' # The problem size is: ', n

  print '(a)', ' #----------------------------------------------------'
  print '(a)', ' # Be careful with memory usage.  If you exceed the'
  print '(a)', ' # 90% (as a rule of thumb) of the memory your GPU'
  print '(a)', ' # card has the program will become memory bound no'
  print '(a)', ' # matter how you wrote it, and it may even crash.'
  print '(a)', ' #----------------------------------------------------'

  print '(a)', ' # Creating a grid'
  call Flow % Create_Field(Grid)

  ! Discretize the matrix for diffusion
  call Process % Form_Diffusion_Matrix(Flow)

  ! Take the alias now
  A => Flow % Nat % M

  allocate(b(n))
  allocate(c(n))

  b(:) = 2.0

  ! Copy operand matrix and vector to the device ...
  ! ... and reserve memory for result vector on device
  call Gpu % Sparse_Copy_To_Device(A)
  call Gpu % Vector_Real_Copy_To_Device(b)
  call Gpu % Vector_Real_Create_On_Device(c)

  !-----------------------------------------------!
  !   Performing a fake time loop on the device   !
  !-----------------------------------------------!
  print '(a,i6,a)', ' # Performing ', N_STEPS, ' sparse-matrix vector products'
  call cpu_time(ts)
  do step = 1, N_STEPS
    call Linalg % Mat_X_Vec(n, c, A, b)
  end do
  call cpu_time(te)

  ! Copy results back to host
  call Gpu % Vector_Update_Host(c)

  ! Destroy data on the device, you don't need them anymore
  call Gpu % Sparse_Destroy_On_Device(A)
  call Gpu % Vector_Real_Destroy_On_Device(b)
  call Gpu % Vector_Real_Destroy_On_Device(c)

  ! Print result
  print '(a,es12.3)', ' vector c(1  ):', c(1)
  print '(a,es12.3)', ' vector c(2  ):', c(2)
  print '(a,es12.3)', ' vector c(n-1):', c(Grid % n_cells - 1)
  print '(a,es12.3)', ' vector c(n  ):', c(Grid % n_cells)

  print '(a,f12.3,a)', ' # Time elapsed for TEST 1: ', te-ts, ' [s]'

  end subroutine

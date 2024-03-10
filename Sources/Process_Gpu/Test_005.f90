!==============================================================================!
  subroutine Test_005
!------------------------------------------------------------------------------!
!>  Tests towards gradient calculation
!------------------------------------------------------------------------------!
  use Field_Mod
  use Gpu_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(Grid_Type)    :: Grid            ! computational grid
  type(Field_Type)   :: Flow            ! flow field
  real, allocatable  :: phi_x(:)        ! gradient in x direction
  real, allocatable  :: phi_y(:)        ! gradient in y direction
  real, allocatable  :: phi_z(:)        ! gradient in z direction
  integer, parameter :: N_STEPS = 120   ! spend enough time on device
  real               :: ts, te
  integer            :: n, c, step
  character(len=11)  :: root_control    = 'control.005'
!==============================================================================!

  print '(a)', ' #===================================================='
  print '(a)', ' # TEST 5: Creating a flow field and gradient matrix'
  print '(a)', ' #===================================================='

  print '(a)', ' # Opening the control file '//root_control
  call Control % Open_Root_File(root_control)

  print '(a)', ' # Creating a grid'
  call Grid % Load_And_Prepare_For_Processing(1)

  n = Grid % n_cells
  print '(a, i12)', ' # The problem size is: ', n

  print '(a)', ' #----------------------------------------------------'
  print '(a)', ' # Be careful with memory usage.  If you exceed the'
  print '(a)', ' # 90% (as a rule of thumb) of the memory your GPU'
  print '(a)', ' # card has the program will become memory bound no'
  print '(a)', ' # matter how you wrote it, and it may even crash.'
  print '(a)', ' #----------------------------------------------------'

  print '(a)', ' # Creating a field'
  call Flow % Create_Field(Grid)

  print '(a)', ' # Allocating and initializing arrays for gradients'
  allocate(phi_x(-Grid % n_bnd_cells:Grid % n_cells));  phi_x(:) = 0.0
  allocate(phi_y(-Grid % n_bnd_cells:Grid % n_cells));  phi_y(:) = 0.0
  allocate(phi_z(-Grid % n_bnd_cells:Grid % n_cells));  phi_z(:) = 0.0

  print '(a)', ' # Initialize phi with something'
  do c = -Grid % n_bnd_cells, Grid % n_cells
    Flow % p % n(c) = 0.111111 * Grid % xc(c)**2  &
                    + 0.222222 * Grid % yc(c)**2  &
                    + 0.333333 * Grid % zc(c)**2
  end do
  call Grid % Save_Debug_Vtu("init", scalar_name="init", scalar_cell=Flow % p % n)

  print '(a)', ' # Calculating gradient matrix for the field'
  call Flow % Calculate_Grad_Matrix()

  ! Copy what you need for gradient calculation to the device
  call Gpu % Matrix_Real_Copy_To_Device(Flow % grad_c2c)
  call Gpu % Vector_Int_Copy_To_Device(Grid % cells_n_cells)
  call Gpu % Matrix_Int_Copy_To_Device(Grid % cells_c)
  call Gpu % Vector_Real_Copy_To_Device(Grid % xc)
  call Gpu % Vector_Real_Copy_To_Device(Grid % yc)
  call Gpu % Vector_Real_Copy_To_Device(Grid % zc)
  call Gpu % Vector_Real_Copy_To_Device(Grid % dx)
  call Gpu % Vector_Real_Copy_To_Device(Grid % dy)
  call Gpu % Vector_Real_Copy_To_Device(Grid % dz)
  call Gpu % Vector_Real_Copy_To_Device(Flow % p % n)
  call Gpu % Vector_Real_Create_On_Device(phi_x)
  call Gpu % Vector_Real_Create_On_Device(phi_y)
  call Gpu % Vector_Real_Create_On_Device(phi_z)

  print '(a,i6,a)', ' # Calculating gradients of the field over ',  &
                    N_STEPS, ' pseudo time steps'
  call cpu_time(ts)
  do step = 1, N_STEPS
    if(mod(step, 12) .eq. 0)  &
      print '(a,i12,es12.3)', ' time step = ', step
    call Flow % Grad_Component(Grid, Flow % p % n, 1, phi_x)
    call Flow % Grad_Component(Grid, Flow % p % n, 2, phi_y)
    call Flow % Grad_Component(Grid, Flow % p % n, 3, phi_z)
  end do
  call cpu_time(te)

  ! Copy results back to host
  call Gpu % Vector_Update_Host(phi_x)
  call Gpu % Vector_Update_Host(phi_y)
  call Gpu % Vector_Update_Host(phi_z)
  call Grid % Save_Debug_Vtu("grad_0", vector_name="graed_0", vector_cell=(/phi_x, phi_y, phi_z/))

  ! Destroy data on the device, you don't need them anymore
  call Gpu % Matrix_Real_Destroy_On_Device(Flow % grad_c2c)
  call Gpu % Vector_Int_Destroy_On_Device(Grid % cells_n_cells)
  call Gpu % Matrix_Int_Destroy_On_Device(Grid % cells_c)
  call Gpu % Vector_Real_Destroy_On_Device(Grid % xc)
  call Gpu % Vector_Real_Destroy_On_Device(Grid % yc)
  call Gpu % Vector_Real_Destroy_On_Device(Grid % zc)
  call Gpu % Vector_Real_Destroy_On_Device(Grid % dx)
  call Gpu % Vector_Real_Destroy_On_Device(Grid % dy)
  call Gpu % Vector_Real_Destroy_On_Device(Grid % dz)
  call Gpu % Vector_Real_Destroy_On_Device(Flow % p % n)
  call Gpu % Vector_Real_Destroy_On_Device(phi_x)
  call Gpu % Vector_Real_Destroy_On_Device(phi_y)
  call Gpu % Vector_Real_Destroy_On_Device(phi_z)

  print '(a,f12.3,a)', ' # Time elapsed for TEST 5: ', te-ts, ' [s]'

  end subroutine

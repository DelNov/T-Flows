#include "../Shared/Browse.h90"

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
  integer, parameter :: N_STEPS = 120   ! spend enough time on device
  integer            :: n, c, step
  character(len=11)  :: root_control = 'control.005'
!==============================================================================!

  ! Start the parallel run and the profiler
  call Global % Start_Parallel
  call Profiler % Start('Test_005')

  O_Print '(a)', ' #===================================================='
  O_Print '(a)', ' # TEST 5: Creating a flow field and gradient matrix'
  O_Print '(a)', ' #===================================================='

  O_Print '(a)', ' # Opening the control file '//root_control
  call Control % Open_Root_File(root_control)

  O_Print '(a)', ' # Creating a grid'
  call Grid % Load_And_Prepare_For_Processing(1)

  n = Grid % n_cells
  O_Print '(a, i12)', ' # The problem size is: ', n

  O_Print '(a)', ' #----------------------------------------------------'
  O_Print '(a)', ' # Be careful with memory usage.  If you exceed the'
  O_Print '(a)', ' # 90% (as a rule of thumb) of the memory your GPU'
  O_Print '(a)', ' # card has the program will become memory bound no'
  O_Print '(a)', ' # matter how you wrote it, and it may even crash.'
  O_Print '(a)', ' #----------------------------------------------------'

  O_Print '(a)', ' # Creating a field'
  call Flow % Create_Field(Grid)

  O_Print '(a)', ' # Initialize phi with something'
  do c = -Grid % n_bnd_cells, Grid % n_cells - Grid % Comm % n_buff_cells
    Flow % p % n(c) = 0.111111 * Grid % xc(c)**2  &
                    + 0.222222 * Grid % yc(c)**2  &
                    + 0.333333 * Grid % zc(c)**2
  end do
  call Grid % Exchange_Cells_Real(Flow % p % n)
  call Grid % Save_Debug_Vtu("init",              &
                             scalar_name="init",  &
                             scalar_cell=Flow % p % n)

  O_Print '(a)', ' # Calculating gradient matrix for the field'
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
  call Gpu % Vector_Real_Create_On_Device(Flow % phi_x)
  call Gpu % Vector_Real_Create_On_Device(Flow % phi_y)
  call Gpu % Vector_Real_Create_On_Device(Flow % phi_z)

  O_Print '(a,i6,a)', ' # Calculating gradients of the field over ',  &
                    N_STEPS, ' pseudo time steps'
  call Profiler % Start('Useful_Work')
  do step = 1, N_STEPS
    if(mod(step, 12) .eq. 0) then
      O_Print '(a,i12,es12.3)', ' time step = ', step
    end if
    call Flow % Grad_Component(Grid, Flow % p % n, 1, Flow % phi_x)
    call Flow % Grad_Component(Grid, Flow % p % n, 2, Flow % phi_y)
    call Flow % Grad_Component(Grid, Flow % p % n, 3, Flow % phi_z)
  end do
  call Profiler % Stop('Useful_Work')

  ! Copy results back to host
  call Gpu % Vector_Update_Host(Flow % phi_x)
  call Gpu % Vector_Update_Host(Flow % phi_y)
  call Gpu % Vector_Update_Host(Flow % phi_z)
  call Grid % Save_Debug_Vtu("grad",                      &
                             vector_name="grad",          &
                             vector_cell=(/Flow % phi_x,  &
                                           Flow % phi_y,  &
                                           Flow % phi_z/))

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
  call Gpu % Vector_Real_Destroy_On_Device(Flow % phi_x)
  call Gpu % Vector_Real_Destroy_On_Device(Flow % phi_y)
  call Gpu % Vector_Real_Destroy_On_Device(Flow % phi_z)

  ! End the profiler and the parallel run
  call Profiler % Stop('Test_005')
  call Profiler % Statistics(indent=24)
  call Global % End_Parallel

  end subroutine

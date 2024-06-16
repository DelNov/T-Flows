!==============================================================================!
  subroutine Grid_Destroy_On_Device(Gpu, Grid, Turb)
!------------------------------------------------------------------------------!
!>  Destroy all the grid variables you don't need in GPU any more.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type) :: Gpu   !! parent class
  type(Grid_Type) :: Grid  !! grid to destroy on device
  type(Turb_Type) :: Turb  !! to check if wall distance should be destroyed
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Gpu)
    Unused(Grid)
    Unused(Turb)
# endif
!==============================================================================!

# if T_FLOWS_GPU == 0
    return
# else
    O_Print '(a)', ' #-------------------------------'
    O_Print '(a)', ' # Destroying grid on the device'
    O_Print '(a)', ' #-------------------------------'
# endif

  call Gpu % Matrix_Int_Destroy_On_Device(Grid % faces_c)
  call Gpu % Vector_Int_Destroy_On_Device(Grid % cells_n_cells)
  call Gpu % Matrix_Int_Destroy_On_Device(Grid % cells_c)
  call Gpu % Matrix_Int_Destroy_On_Device(Grid % cells_f)
  call Gpu % Vector_Real_Destroy_On_Device(Grid % xc)
  call Gpu % Vector_Real_Destroy_On_Device(Grid % yc)
  call Gpu % Vector_Real_Destroy_On_Device(Grid % zc)
  call Gpu % Vector_Real_Destroy_On_Device(Grid % dx)
  call Gpu % Vector_Real_Destroy_On_Device(Grid % dy)
  call Gpu % Vector_Real_Destroy_On_Device(Grid % dz)
  call Gpu % Vector_Real_Destroy_On_Device(Grid % sx)
  call Gpu % Vector_Real_Destroy_On_Device(Grid % sy)
  call Gpu % Vector_Real_Destroy_On_Device(Grid % sz)
  call Gpu % Vector_Real_Destroy_On_Device(Grid % s)
  call Gpu % Vector_Real_Destroy_On_Device(Grid % d)
  call Gpu % Vector_Real_Destroy_On_Device(Grid % vol)
  call Gpu % Vector_Int_Destroy_On_Device(Grid % region % f_face)
  call Gpu % Vector_Int_Destroy_On_Device(Grid % region % l_face)
  call Gpu % Vector_Int_Destroy_On_Device(Grid % region % f_cell)
  call Gpu % Vector_Int_Destroy_On_Device(Grid % region % l_cell)
  if(Turb % model .ne. NO_TURBULENCE_MODEL) then
    call Gpu % Vector_Real_Destroy_On_Device(Grid % wall_dist)
  end if

  end subroutine


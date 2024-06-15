!==============================================================================!
  subroutine Grid_Copy_To_Device(Gpu, Grid)
!------------------------------------------------------------------------------!
!>  Copy all the grid variables you need in your simulation to GPU.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type) :: Gpu   !! parent class
  type(Grid_Type) :: Grid  !! grid to transfer to device
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Gpu)
    Unused(Grid)
# endif
!==============================================================================!

# if T_FLOWS_GPU == 0
    return
# else
    O_Print '(a)', ' #----------------------------'
    O_Print '(a)', ' # Copying grid to the device'
    O_Print '(a)', ' #----------------------------'
# endif

  call Gpu % Matrix_Int_Copy_To_Device(Grid % faces_c)
  call Gpu % Vector_Int_Copy_To_Device(Grid % cells_n_cells)
  call Gpu % Matrix_Int_Copy_To_Device(Grid % cells_c)
  call Gpu % Matrix_Int_Copy_To_Device(Grid % cells_f)
  call Gpu % Vector_Real_Copy_To_Device(Grid % xc)
  call Gpu % Vector_Real_Copy_To_Device(Grid % yc)
  call Gpu % Vector_Real_Copy_To_Device(Grid % zc)
  call Gpu % Vector_Real_Copy_To_Device(Grid % dx)
  call Gpu % Vector_Real_Copy_To_Device(Grid % dy)
  call Gpu % Vector_Real_Copy_To_Device(Grid % dz)
  call Gpu % Vector_Real_Copy_To_Device(Grid % sx)
  call Gpu % Vector_Real_Copy_To_Device(Grid % sy)
  call Gpu % Vector_Real_Copy_To_Device(Grid % sz)
  call Gpu % Vector_Real_Copy_To_Device(Grid % s)
  call Gpu % Vector_Real_Copy_To_Device(Grid % d)
  call Gpu % Vector_Real_Copy_To_Device(Grid % vol)
  call Gpu % Vector_Int_Copy_To_Device(Grid % region % f_face)
  call Gpu % Vector_Int_Copy_To_Device(Grid % region % l_face)
  call Gpu % Vector_Int_Copy_To_Device(Grid % region % f_cell)
  call Gpu % Vector_Int_Copy_To_Device(Grid % region % l_cell)

  end subroutine


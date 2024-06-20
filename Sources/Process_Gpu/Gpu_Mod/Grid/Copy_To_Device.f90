!==============================================================================!
  subroutine Grid_Copy_To_Device(Gpu, Turb, Grid)
!------------------------------------------------------------------------------!
!>  Copy all the grid variables you need in your simulation to GPU.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type)         :: Gpu   !! parent class
  type(Turb_Type)         :: Turb  !! to check if wall distance is needed
  type(Grid_Type), target :: Grid  !! grid to transfer to device
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Gpu)
    Unused(Turb)
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
  grid_faces_c       => Grid % faces_c
  grid_cells_n_cells => Grid % cells_n_cells
  grid_cells_c       => Grid % cells_c
  grid_cells_f       => Grid % cells_f
  call Gpu % Vector_Real_Copy_To_Device(Grid % xc)
  call Gpu % Vector_Real_Copy_To_Device(Grid % yc)
  call Gpu % Vector_Real_Copy_To_Device(Grid % zc)
  grid_xc  => Grid % xc
  grid_yc  => Grid % yc
  grid_zc  => Grid % zc
  call Gpu % Vector_Real_Copy_To_Device(Grid % dx)
  call Gpu % Vector_Real_Copy_To_Device(Grid % dy)
  call Gpu % Vector_Real_Copy_To_Device(Grid % dz)
  grid_dx  => Grid % dx
  grid_dy  => Grid % dy
  grid_dz  => Grid % dz
  call Gpu % Vector_Real_Copy_To_Device(Grid % sx)
  call Gpu % Vector_Real_Copy_To_Device(Grid % sy)
  call Gpu % Vector_Real_Copy_To_Device(Grid % sz)
  grid_sx  => Grid % sx
  grid_sy  => Grid % sy
  grid_sz  => Grid % sz
  call Gpu % Vector_Real_Copy_To_Device(Grid % s)
  call Gpu % Vector_Real_Copy_To_Device(Grid % d)
  grid_s   => Grid % s
  grid_d   => Grid % d
  call Gpu % Vector_Real_Copy_To_Device(Grid % vol)
  grid_vol => Grid % vol
  call Gpu % Vector_Int_Copy_To_Device(Grid % region % f_face)
  call Gpu % Vector_Int_Copy_To_Device(Grid % region % l_face)
  call Gpu % Vector_Int_Copy_To_Device(Grid % region % f_cell)
  call Gpu % Vector_Int_Copy_To_Device(Grid % region % l_cell)
  grid_n_regions     =  Grid % n_regions
  grid_n_bnd_regions =  Grid % n_bnd_regions
  grid_region_f_face => Grid % region % f_face
  grid_region_l_face => Grid % region % l_face
  grid_region_f_cell => Grid % region % f_cell
  grid_region_l_cell => Grid % region % l_cell
  if(Turb % model .ne. NO_TURBULENCE_MODEL) then
    call Gpu % Vector_Real_Copy_To_Device(Grid % wall_dist)
    grid_wall_dist => Grid % wall_dist
  end if

  end subroutine


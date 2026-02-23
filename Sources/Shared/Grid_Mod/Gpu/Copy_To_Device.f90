!==============================================================================!
  subroutine Copy_Grid_To_Device(Grid)
!------------------------------------------------------------------------------!
!>  Copy all the grid variables you need in your simulation to GPU.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type), target :: Grid  !! parent grid to transfer to device
!-----------------------------------[Locals]-----------------------------------!
  integer :: lb, ub
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Grid)
# endif
!==============================================================================!

# if T_FLOWS_GPU == 1
    O_Print '(a)', ' #----------------------------'
    O_Print '(a)', ' # Copying grid to the device'
    O_Print '(a)', ' #----------------------------'
# endif

  call Gpu % Matrix_Int_Copy_To_Device(Grid % faces_c)
  lb = lbound(Grid % cells_n_cells,1)
  ub = ubound(Grid % cells_n_cells,1)
  call Gpu % Vector_Int_Copy_To_Device(lb,ub, Grid % cells_n_cells)
  lb = lbound(Grid % cells_i_cells,1)
  ub = ubound(Grid % cells_i_cells,1)
  call Gpu % Vector_Int_Copy_To_Device(lb,ub, Grid % cells_i_cells)
  call Gpu % Matrix_Int_Copy_To_Device(Grid % cells_c)
  call Gpu % Matrix_Int_Copy_To_Device(Grid % cells_f)
  grid_faces_c       => Grid % faces_c
  grid_cells_n_cells => Grid % cells_n_cells
  grid_cells_i_cells => Grid % cells_i_cells
  grid_cells_c       => Grid % cells_c
  grid_cells_f       => Grid % cells_f
  lb = lbound(Grid % xc,1)
  ub = ubound(Grid % xc,1)
  call Gpu % Vector_Real_Copy_To_Device(lb,ub, Grid % xc)
  lb = lbound(Grid % yc,1)
  ub = ubound(Grid % yc,1)
  call Gpu % Vector_Real_Copy_To_Device(lb,ub, Grid % yc)
  lb = lbound(Grid % zc,1)
  ub = ubound(Grid % zc,1)
  call Gpu % Vector_Real_Copy_To_Device(lb,ub, Grid % zc)
  grid_xc  => Grid % xc
  grid_yc  => Grid % yc
  grid_zc  => Grid % zc
  lb = lbound(Grid % dx,1)
  ub = ubound(Grid % dx,1)
  call Gpu % Vector_Real_Copy_To_Device(lb,ub, Grid % dx)
  lb = lbound(Grid % dy,1)
  ub = ubound(Grid % dy,1)
  call Gpu % Vector_Real_Copy_To_Device(lb,ub, Grid % dy)
  lb = lbound(Grid % dz,1)
  ub = ubound(Grid % dz,1)
  call Gpu % Vector_Real_Copy_To_Device(lb,ub, Grid % dz)
  grid_dx  => Grid % dx
  grid_dy  => Grid % dy
  grid_dz  => Grid % dz
  lb = lbound(Grid % sx,1)
  ub = ubound(Grid % sx,1)
  call Gpu % Vector_Real_Copy_To_Device(lb,ub, Grid % sx)
  lb = lbound(Grid % sy,1)
  ub = ubound(Grid % sy,1)
  call Gpu % Vector_Real_Copy_To_Device(lb,ub, Grid % sy)
  lb = lbound(Grid % sz,1)
  ub = ubound(Grid % sz,1)
  call Gpu % Vector_Real_Copy_To_Device(lb,ub, Grid % sz)
  grid_sx  => Grid % sx
  grid_sy  => Grid % sy
  grid_sz  => Grid % sz
  lb = lbound(Grid % f,1)
  ub = ubound(Grid % f,1)
  call Gpu % Vector_Real_Copy_To_Device(lb,ub, Grid % f)
  lb = lbound(Grid % s,1)
  ub = ubound(Grid % s,1)
  call Gpu % Vector_Real_Copy_To_Device(lb,ub, Grid % s)
  lb = lbound(Grid % d,1)
  ub = ubound(Grid % d,1)
  call Gpu % Vector_Real_Copy_To_Device(lb,ub, Grid % d)
  grid_f   => Grid % f
  grid_s   => Grid % s
  grid_d   => Grid % d
  lb = lbound(Grid % vol,1)
  ub = ubound(Grid % vol,1)
  call Gpu % Vector_Real_Copy_To_Device(lb,ub, Grid % vol)
  grid_vol => Grid % vol
  lb = lbound(Grid % region % f_face,1)
  ub = ubound(Grid % region % f_face,1)
  call Gpu % Vector_Int_Copy_To_Device(lb,ub, Grid % region % f_face)
  lb = lbound(Grid % region % l_face,1)
  ub = ubound(Grid % region % l_face,1)
  call Gpu % Vector_Int_Copy_To_Device(lb,ub, Grid % region % l_face)
  lb = lbound(Grid % region % f_cell,1)
  ub = ubound(Grid % region % f_cell,1)
  call Gpu % Vector_Int_Copy_To_Device(lb,ub, Grid % region % f_cell)
  lb = lbound(Grid % region % l_cell,1)
  ub = ubound(Grid % region % l_cell,1)
  call Gpu % Vector_Int_Copy_To_Device(lb,ub, Grid % region % l_cell)
  grid_n_regions     =  Grid % n_regions
  grid_n_bnd_regions =  Grid % n_bnd_regions
  grid_region_f_face => Grid % region % f_face
  grid_region_l_face => Grid % region % l_face
  grid_region_f_cell => Grid % region % f_cell
  grid_region_l_cell => Grid % region % l_cell
  lb = lbound(Grid % wall_dist,1)
  ub = ubound(Grid % wall_dist,1)
  call Gpu % Vector_Real_Copy_To_Device(lb,ub, Grid % wall_dist)
  grid_wall_dist => Grid % wall_dist

  end subroutine


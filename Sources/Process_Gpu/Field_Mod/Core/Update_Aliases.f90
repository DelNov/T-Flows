!==============================================================================!
  subroutine Update_Aliases(Flow)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow  !! parent flow object
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
!==============================================================================!

  ! Aliases connected to Flow
  u_n => Flow % u % n
  v_n => Flow % v % n
  w_n => Flow % w % n

  u_o => Flow % u % o
  v_o => Flow % v % o
  w_o => Flow % w % o

  pp_n => Flow % pp % n
  pp_x => Flow % pp % x
  pp_y => Flow % pp % y
  pp_z => Flow % pp % z

  p_n => Flow % p % n
  p_x => Flow % p % x
  p_y => Flow % p % y
  p_z => Flow % p % z

  ! Aliases connected to Grid
  Grid => Flow % pnt_grid

  grid_n_cells      = Grid % n_cells
  grid_n_faces      = Grid % n_faces
  grid_n_bnd_cells  = Grid % n_bnd_cells
  grid_n_buff_cells = Grid % Comm % n_buff_cells
  grid_n_regions    = Grid % n_regions

  grid_faces_c       => Grid % faces_c
  grid_cells_c       => Grid % cells_c
  grid_cells_f       => Grid % cells_f
  grid_cells_n_cells => Grid % cells_n_cells

  grid_reg_f_face => Grid % region % f_face
  grid_reg_l_face => Grid % region % l_face

  grid_vol => Grid % vol
  grid_xc  => Grid % xc
  grid_yc  => Grid % yc
  grid_zc  => Grid % zc
  grid_dx  => Grid % dx
  grid_dy  => Grid % dy
  grid_dz  => Grid % dz
  grid_sx  => Grid % sx
  grid_sy  => Grid % sy
  grid_sz  => Grid % sz

  end subroutine

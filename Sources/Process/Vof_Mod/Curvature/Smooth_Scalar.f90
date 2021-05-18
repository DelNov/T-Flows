!==============================================================================!
  subroutine Smooth_Scalar(Vof, Grid, var, smooth_var, n_conv)
!------------------------------------------------------------------------------!
!   Smoothes scalar using a Laplacian smoother                                 !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: sum_vol_area => r_cell_01,   &
                      sum_area     => r_cell_02
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  type(Grid_Type)         :: Grid
  integer                 :: n_conv, nb, nc
  real                    :: var(-Grid % n_bnd_cells: Grid % n_cells)
  real                    :: smooth_var(-Grid % n_bnd_cells    &
                                       : Grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, c, c1, c2, c_iter
  real    :: fs, vol_face
!==============================================================================!

  ! Take aliases
  nb = Grid % n_bnd_cells
  nc = Grid % n_cells

  ! Copy the values from phi % n to local variable
  smooth_var(:) = var(:)

  do c_iter = 1, n_conv

    sum_vol_area(-nb:nc) = 0.0
    sum_area    (-nb:nc) = 0.0

    !-------------------------------!
    !   Extrapolate to boundaries   !
    !-------------------------------!

    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      if(c2 < 0) smooth_var(c2) = smooth_var(c1)
    end do

    ! At boundaries
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      if(c2 < 0) then
        sum_vol_area(c1) = sum_vol_area(c1) + smooth_var(c1) * Grid % s(s)
        sum_area(c1) = sum_area(c1) + Grid % s(s)
      end if
    end do

    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      if(c2 > 0) then
        fs = Grid % f(s)
        vol_face = fs * smooth_var(c1) + (1.0 - fs) * smooth_var(c2)
        sum_vol_area(c1) = sum_vol_area(c1) + vol_face * Grid % s(s)
        sum_vol_area(c2) = sum_vol_area(c2) + vol_face * Grid % s(s)
        sum_area(c1) = sum_area(c1) + Grid % s(s)
        sum_area(c2) = sum_area(c2) + Grid % s(s)
      end if
    end do

    call Grid % Exchange_Cells_Real(sum_vol_area(-nb:nc))
    call Grid % Exchange_Cells_Real(sum_area    (-nb:nc))

    do c = 1, Grid % n_cells
      smooth_var(c) = max(min(sum_vol_area(c) / sum_area(c), 1.0), 0.0)
      smooth_var(c) = sum_vol_area(c) / sum_area(c)
    end do
    call Grid % Exchange_Cells_Real(smooth_var)

  end do

  ! At boundaries
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2 < 0) smooth_var(c2) = smooth_var(c1)
  end do

  call Grid % Exchange_Cells_Real(smooth_var)

  end subroutine

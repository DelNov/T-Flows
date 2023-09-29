!==============================================================================!
  subroutine Smooth_Scalar(Vof, Grid, var, smooth_var, n_conv)
!------------------------------------------------------------------------------!
!   Smoothes scalar using a Laplacian smoother                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  type(Grid_Type)         :: Grid
  integer                 :: n_conv, nb, nc, reg
  real                    :: var(-Grid % n_bnd_cells: Grid % n_cells)
  real                    :: smooth_var(-Grid % n_bnd_cells    &
                                       : Grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer                   :: s, c, c1, c2, c_iter
  real                      :: fs, vol_face
  real, contiguous, pointer :: sum_vol_area(:), sum_area(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Vof)
!==============================================================================!

  call Work % Connect_Real_Cell(sum_vol_area, sum_area)

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
    do reg = Boundary_Regions()
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)
        smooth_var(c2) = smooth_var(c1)
      end do
    end do

    ! At boundaries
    do reg = Boundary_Regions()
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        sum_vol_area(c1) = sum_vol_area(c1) + smooth_var(c1) * Grid % s(s)
        sum_area(c1) = sum_area(c1) + Grid % s(s)
      end do
    end do

    do s = Faces_In_Domain_And_At_Buffers()
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      fs = Grid % f(s)
      vol_face         = fs * smooth_var(c1) + (1.0-fs) * smooth_var(c2)
      sum_vol_area(c1) = sum_vol_area(c1) + vol_face * Grid % s(s)
      sum_vol_area(c2) = sum_vol_area(c2) + vol_face * Grid % s(s)
      sum_area(c1)     = sum_area(c1) + Grid % s(s)
      sum_area(c2)     = sum_area(c2) + Grid % s(s)
    end do

    call Grid % Exchange_Cells_Real(sum_vol_area(-nb:nc))
    call Grid % Exchange_Cells_Real(sum_area    (-nb:nc))

    do c = Cells_In_Domain()
      smooth_var(c) = max(min(sum_vol_area(c) / sum_area(c), 1.0), 0.0)
      smooth_var(c) = sum_vol_area(c) / sum_area(c)
    end do
    call Grid % Exchange_Cells_Real(smooth_var)

  end do

  ! At boundaries
  do reg = Boundary_Regions()
    do s = Faces_In_Region(reg)
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      smooth_var(c2) = smooth_var(c1)
    end do
  end do
  call Grid % Exchange_Cells_Real(smooth_var(-nb:nc))

  call Work % Disconnect_Real_Cell(sum_vol_area, sum_area)

  end subroutine

!==============================================================================!
  subroutine Multiphase_Mod_Vof_Smooth_Scalar(grid, mult, var,   &
                                              smooth_var, n_conv)
!------------------------------------------------------------------------------!
!    Smoothes scalar using a laplacian smoother                                !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: sum_vol_area => r_cell_01,   &
                      sum_area     => r_cell_02
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Grid_Type)               :: grid
  integer                       :: n_conv, nb, nc
  real                          :: var(-grid % n_bnd_cells: grid % n_cells)
  real                          :: smooth_var(-grid % n_bnd_cells    &
                                             : grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type), pointer :: vof
  integer                 :: s, c, c1, c2, c_iter
  real                    :: fs, vol_face
!==============================================================================!

  ! Take aliases
  vof => mult % vof

  nb = grid % n_bnd_cells
  nc = grid % n_cells

  ! Copy the values from phi % n to local variable
  smooth_var(:) = var(:)

  do c_iter = 1, n_conv

    sum_vol_area(-nb:nc) = 0.0
    sum_area    (-nb:nc) = 0.0

    !-------------------------------!
    !   Extrapolate to boundaries   !
    !-------------------------------!

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 < 0) smooth_var(c2) = smooth_var(c1)
    end do

    ! At boundaries
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 < 0) then
        sum_vol_area(c1) = sum_vol_area(c1) + smooth_var(c1) * grid % s(s)
        sum_area(c1) = sum_area(c1) + grid % s(s)
      end if
    end do

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      if(c2 > 0) then
        fs = grid % f(s)
        vol_face = fs * smooth_var(c1) + (1.0 - fs) * smooth_var(c2)
        sum_vol_area(c1) = sum_vol_area(c1) + vol_face * grid % s(s)
        sum_vol_area(c2) = sum_vol_area(c2) + vol_face * grid % s(s)
        sum_area(c1) = sum_area(c1) + grid % s(s)
        sum_area(c2) = sum_area(c2) + grid % s(s)
      end if
    end do

    call Grid_Mod_Exchange_Cells_Real(grid, sum_vol_area(-nb:nc))
    call Grid_Mod_Exchange_Cells_Real(grid, sum_area    (-nb:nc))

    do c = 1, grid % n_cells
      smooth_var(c) = max(min(sum_vol_area(c) / sum_area(c), 1.0), 0.0)
      smooth_var(c) = sum_vol_area(c) / sum_area(c)
    end do
    call Grid_Mod_Exchange_Cells_Real(grid, smooth_var)

  end do

  ! At boundaries
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 < 0) smooth_var(c2) = smooth_var(c1)
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, smooth_var)

  end subroutine

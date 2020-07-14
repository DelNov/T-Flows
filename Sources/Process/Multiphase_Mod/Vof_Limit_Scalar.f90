!==============================================================================!
  subroutine Multiphase_Mod_Vof_Limit_Scalar(grid, mult, var, limit_var)
!------------------------------------------------------------------------------!
!   Smoothes scalar using a laplacian smoother                                 !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: w             => r_cell_01,   &
                      w_int         => r_cell_02,   &
                      face_scalar   => r_face_01,   &
                      wf            => r_face_02
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Grid_Type)               :: grid
  integer                       :: n_conv
  real                          :: var(-grid % n_bnd_cells:grid % n_cells)
  real                          :: limit_var(-grid % n_bnd_cells:grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type),    pointer :: vof
  type(Matrix_Type), pointer :: a
  integer                    :: i, s, c, c1, c2, c_iter, run
  integer                    :: face_init, face_end, face_step
  real                       :: fs, vol_face, csk, w1, w2, epsloc
!==============================================================================!

  ! Take aliases
  if (mult % d_func) then
    vof => mult % dist_func
  else
    vof => mult % vof
  end if

  epsloc = epsilon(epsloc)

  ! Copy the values from phi % n to local variable
  limit_var(:) = var(:)

  csk = 0.5
  i = 1

  ! find w

  do c = 1, grid % n_cells
    w(c) = sqrt(vof % n(c) * (1.0 - vof % n(c)) + MICRO)
  end do

  do c_iter = 0, i

    !-------------------------------!
    !   Extrapolate to boundaries   !
    !-------------------------------!

    ! At boundaries
    do s = 1, grid % n_bnd_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      limit_var(c2) = limit_var(c1)
      face_scalar(s) = limit_var(c2) * w(c1)
      wf(s) = w(c1)
    end do

    ! Interior faces
    do s = grid % n_bnd_faces + 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)
      face_scalar(s) = fs * limit_var(c1) * w(c1)           &
                     + (1.0 - fs) * limit_var(c2) * w(c2)
      wf(s) = fs * w(c1) + (1.0 - fs) * w(c2)
    end do

    ! Back to cells
    call Multiphase_Mod_Vof_Interpolate_Faces_Cells(grid,                     &
                                                    face_scalar, limit_var)

    call Multiphase_Mod_Vof_Interpolate_Faces_Cells(grid, wf, w_int)

    do c = 1, grid % n_cells
      limit_var(c) = 2.0 * sqrt(vof % n(c) * (1.0 - vof % n(c))) * var(c)   &
                   + (1.0 - 2.0 * sqrt(vof % n(c) * (1.0 - vof % n(c))))    &
                   * limit_var(c) / w_int(c)
    end do

    call Grid_Mod_Exchange_Cells_Real(grid, limit_var)
  end do

  ! At boundaries
  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    limit_var(c2) = limit_var(c1)
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, limit_var)

  end subroutine

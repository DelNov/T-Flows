!==============================================================================!
  subroutine Multiphase_Mod_Vof_Find_Weight_Faces(grid)
!------------------------------------------------------------------------------!
!   Computes weights for cell interpolation from faces                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, s, i_fac   ! counters
  real    :: rx, ry, rz, lambda_x, lambda_y, lambda_z
  real    :: ixx, iyy, izz, ixz, iyz, ixy, d
  real    :: a11, a12, a13, a21, a22, a23, a31, a32, a33
  real    :: corr_x, corr_y, corr_z
  real    :: xf, yf, zf
!==============================================================================!

  !---------------------------------------!
  !   Browse through all faces to form    !
  !              weight_faces             !
  !---------------------------------------!

  ! first find the limits of boundary faces
  ! loops are nicer this way
  loop_b_faces: do s = 1, grid % n_faces
    c2 = grid % faces_c(2,s)
    if (c2 > 0) then
      grid % n_bnd_faces = s - 1
      exit loop_b_faces
    end if
  end do loop_b_faces

  do c = 1, grid % n_cells
    rx = 0.0; ry = 0.0; rz = 0.0
    ixx = 0.0; iyy = 0.0; izz = 0.0; ixz = 0.0; iyz = 0.0; ixy = 0.0
    ! Loop on faces
    do i_fac = 1, grid % cells_n_faces(c)
      s = grid % cells_f(i_fac, c)
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c == c1) then
        rx = rx + (grid % xf(s) - grid % xc(c))
        ry = ry + (grid % yf(s) - grid % yc(c))
        rz = rz + (grid % zf(s) - grid % zc(c))
        ixx = ixx + (grid % xf(s) - grid % xc(c)) ** 2
        iyy = iyy + (grid % yf(s) - grid % yc(c)) ** 2
        izz = izz + (grid % zf(s) - grid % zc(c)) ** 2
        ixy = ixy + (grid % xf(s) - grid % xc(c)) * (grid % yf(s) - grid % yc(c))
        ixz = ixz + (grid % xf(s) - grid % xc(c)) * (grid % zf(s) - grid % zc(c))
        iyz = iyz + (grid % yf(s) - grid % yc(c)) * (grid % zf(s) - grid % zc(c))
      else

        ! Correction for periodic faces:
        call Grid_Mod_Correction_Periodicity(grid, s,   &
                                             corr_x, corr_y, corr_z)
        xf = grid % xf(s) + corr_x
        yf = grid % yf(s) + corr_y
        zf = grid % zf(s) + corr_z

        rx = rx + (xf - grid % xc(c))
        ry = ry + (yf - grid % yc(c))
        rz = rz + (zf - grid % zc(c))
        ixx = ixx + (xf - grid % xc(c)) ** 2
        iyy = iyy + (yf - grid % yc(c)) ** 2
        izz = izz + (zf - grid % zc(c)) ** 2
        ixy = ixy + (xf - grid % xc(c)) * (yf - grid % yc(c))
        ixz = ixz + (xf - grid % xc(c)) * (zf - grid % zc(c))
        iyz = iyz + (yf - grid % yc(c)) * (zf - grid % zc(c))
      end if
    end do
    a11 = iyy * izz - iyz ** 2
    a12 = ixz * iyz - ixy * izz
    a13 = ixy * iyz - ixz * iyy
    a21 = ixz * iyz - ixy * izz
    a22 = ixx * izz - ixz ** 2
    a23 = ixy * ixz - ixx * iyz
    a31 = ixy * iyz - ixz * iyy
    a32 = ixz * ixy - ixx * iyz
    a33 = ixx * iyy - ixy ** 2

    d = ixx * iyy * izz - ixx * iyz **2 - iyy * ixz ** 2 - izz * ixy ** 2   &
      + 2.0 * ixy * ixz * iyz

    lambda_x = (rx * a11 + ry * a12 + rz * a13) / d
    lambda_y = (rx * a21 + ry * a22 + rz * a23) / d
    lambda_z = (rx * a31 + ry * a32 + rz * a33) / d

    do i_fac = 1, grid % cells_n_faces(c)
      s = grid % cells_f(i_fac, c)
      grid % weight_faces(i_fac, c) = 1.0                                 &
                             + lambda_x * (grid % xf(s) - grid % xc(c))   &
                             + lambda_y * (grid % yf(s) - grid % yc(c))   &
                             + lambda_z * (grid % zf(s) - grid % zc(c))
    end do
  end do

  end subroutine

!==============================================================================!
  subroutine Multiphase_Mod_Vof_Find_Weight_Cells_To_Nodes(grid)
!------------------------------------------------------------------------------!
!   Computes weights for cell interpolation to nodes                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer         :: c, n, i_cell, c1, c2, s, i_fac   ! counters
  real            :: rx, ry, rz, lambda_x, lambda_y, lambda_z
  real            :: ixx, iyy, izz, ixz, iyz, ixy, d
  real            :: a11, a12, a13, a21, a22, a23, a31, a32, a33
  real            :: corr_x, corr_y, corr_z
  real            :: epsloc
!==============================================================================!

  epsloc = epsilon(epsloc)

  !---------------------------------------!
  !   Browse through all nodes to form    !
  !              node weights             !
  !---------------------------------------!
  allocate(grid % nodes_weight_c(size(grid % nodes_c,1),size(grid % nodes_c,2)))
  grid % nodes_weight_c = 0.0

  ! Loop on nodes
  do n = 1, grid % n_nodes
    rx = 0.0; ry = 0.0; rz = 0.0
    ixx = 0.0; iyy = 0.0; izz = 0.0; ixz = 0.0; iyz = 0.0; ixy = 0.0
    ! Loop on cells
    do i_cell = 1, grid % nodes_n_cells(n)
      c = grid % nodes_c(i_cell, n)
      rx = rx + (grid % xc(c) - grid % xn(n))
      ry = ry + (grid % yc(c) - grid % yn(n))
      rz = rz + (grid % zc(c) - grid % zn(n))
      ixx = ixx + (grid % xc(c) - grid % xn(n)) ** 2
      iyy = iyy + (grid % yc(c) - grid % yn(n)) ** 2
      izz = izz + (grid % zc(c) - grid % zn(n)) ** 2
      ixy = ixy + (grid % xc(c) - grid % xn(n)) * (grid % yc(c) - grid % yn(n))
      ixz = ixz + (grid % xc(c) - grid % xn(n)) * (grid % zc(c) - grid % zn(n))
      iyz = iyz + (grid % yc(c) - grid % yn(n)) * (grid % zc(c) - grid % zn(n))
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

    lambda_x = (rx * a11 + ry * a12 + rz * a13) / (d + epsloc)
    lambda_y = (rx * a21 + ry * a22 + rz * a23) / (d + epsloc)
    lambda_z = (rx * a31 + ry * a32 + rz * a33) / (d + epsloc)

    do i_cell = 1, grid % nodes_n_cells(n)
      c = grid % nodes_c(i_cell, n)
      grid % nodes_weight_c(i_cell, n) = 1.0                               &
                              + lambda_x * (grid % xc(c) - grid % xn(n))   &
                              + lambda_y * (grid % yc(c) - grid % yn(n))   &
                              + lambda_z * (grid % zc(c) - grid % zn(n))

    end do

  end do

  end subroutine

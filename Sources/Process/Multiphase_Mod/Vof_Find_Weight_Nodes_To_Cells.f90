!==============================================================================!
  subroutine Multiphase_Mod_Vof_Find_Weight_Nodes_To_Cells(grid)
!------------------------------------------------------------------------------!
!   Computes weights for node interpolation to cells                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer         :: c, n, i_nod, c1, c2, s, i_fac   ! counters
  real            :: rx, ry, rz, lambda_x, lambda_y, lambda_z
  real            :: ixx, iyy, izz, ixz, iyz, ixy, d
  real            :: a11, a12, a13, a21, a22, a23, a31, a32, a33
  real            :: corr_x, corr_y, corr_z
  real            :: xf, yf, zf
  real            :: epsloc
!==============================================================================!

  epsloc = epsilon(epsloc)

  !---------------------------------------!
  !   Browse through all nodes to form    !
  !              node weights             !
  !---------------------------------------!
  allocate(grid % cells_weight_n(size(grid % cells_n,1),size(grid % cells_n,2)))
  grid % cells_weight_n= 0.0

  ! Loop on cells
  do c = 1, grid % n_cells
    rx = 0.0; ry = 0.0; rz = 0.0
    ixx = 0.0; iyy = 0.0; izz = 0.0; ixz = 0.0; iyz = 0.0; ixy = 0.0
    ! Loop on nodes
    do i_nod = 1, grid % cells_n_nodes(c)
      n = grid % cells_n(i_nod, c)
      rx = rx + (grid % xn(n) - grid % xc(c))
      ry = ry + (grid % yn(n) - grid % yc(c))
      rz = rz + (grid % zn(n) - grid % zc(c))
      ixx = ixx + (grid % xn(n) - grid % xc(c)) ** 2
      iyy = iyy + (grid % yn(n) - grid % yc(c)) ** 2
      izz = izz + (grid % zn(n) - grid % zc(c)) ** 2
      ixy = ixy + (grid % xn(n) - grid % xc(c)) * (grid % yn(n) - grid % yc(c))
      ixz = ixz + (grid % xn(n) - grid % xc(c)) * (grid % zn(n) - grid % zc(c))
      iyz = iyz + (grid % yn(n) - grid % yc(c)) * (grid % zn(n) - grid % zc(c))
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

    do i_nod = 1, grid % cells_n_nodes(c)
      n = grid % cells_n(i_nod, c)
      grid % cells_weight_n(i_nod, c) = 1.0                               &
                             + lambda_x * (grid % xn(n) - grid % xc(c))   &
                             + lambda_y * (grid % yn(n) - grid % yc(c))   &
                             + lambda_z * (grid % zn(n) - grid % zc(c))
    end do
  end do

  end subroutine

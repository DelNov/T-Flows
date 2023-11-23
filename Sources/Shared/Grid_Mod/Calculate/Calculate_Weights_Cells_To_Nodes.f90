!==============================================================================!
  subroutine Calculate_Weights_Cells_To_Nodes(Grid)
!------------------------------------------------------------------------------!
!>  Computes weights for interpolation from cell to nodes.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid  !! grid under consideration
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, n, i_cell
  real    :: lx, ly, lz, rx, ry, rz, lambda_x, lambda_y, lambda_z
  real    :: ixx, iyy, izz, ixz, iyz, ixy, d
  real    :: a11, a12, a13, a21, a22, a23, a31, a32, a33
  real    :: tot
!==============================================================================!

  ! Allocate memory
  allocate(Grid % weight_c2n(size(Grid % nodes_c,1),size(Grid % nodes_c,2)))
  Grid % weight_c2n = 0.0

  !-------------------------------------------!
  !   Browse through all nodes to calculate   !
  !     weights from cells surrounding it     !
  !-------------------------------------------!
  do n = 1, Grid % n_nodes
    rx = 0.0; ry = 0.0; rz = 0.0
    ixx = 0.0; iyy = 0.0; izz = 0.0; ixz = 0.0; iyz = 0.0; ixy = 0.0

    ! Loop on cells surrounding the node
    do i_cell = 1, Grid % nodes_n_cells(n)
      c = Grid % nodes_c(i_cell, n)

      ! Default distances
      lx = Grid % xc(c) - Grid % xn(n)
      ly = Grid % yc(c) - Grid % yn(n)
      lz = Grid % zc(c) - Grid % zn(n)

      ! Check if periodicity in all three directions and correct distances
      if( abs(lx - Grid % per_x) < abs(lx) ) then
        lx = lx - Grid % per_x
      else if( abs(lx + Grid % per_x) < abs(lx) ) then
        lx = lx + Grid % per_x
      end if
      if( abs(ly - Grid % per_y) < abs(ly) ) then
        ly = ly - Grid % per_y
      else if( abs(ly + Grid % per_y) < abs(ly) ) then
        ly = ly + Grid % per_y
      end if
      if( abs(lz - Grid % per_z) < abs(lz) ) then
        lz = lz - Grid % per_z
      else if( abs(lz + Grid % per_z) < abs(lz) ) then
        lz = lz + Grid % per_z
      end if

      rx  = rx  + lx
      ry  = ry  + ly
      rz  = rz  + lz
      ixx = ixx + lx ** 2
      iyy = iyy + ly ** 2
      izz = izz + lz ** 2
      ixy = ixy + lx * ly
      ixz = ixz + lx * lz
      iyz = iyz + ly * lz
    end do  ! through cells surrounding the node

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

    lambda_x = (rx * a11 + ry * a12 + rz * a13) / (d + FEMTO)
    lambda_y = (rx * a21 + ry * a22 + rz * a23) / (d + FEMTO)
    lambda_z = (rx * a31 + ry * a32 + rz * a33) / (d + FEMTO)

    do i_cell = 1, Grid % nodes_n_cells(n)
      c  = Grid % nodes_c(i_cell, n)

      ! Default distances
      lx = Grid % xc(c) - Grid % xn(n)
      ly = Grid % yc(c) - Grid % yn(n)
      lz = Grid % zc(c) - Grid % zn(n)

      ! Check if periodicity in all three directions and correct distances
      if( abs(lx - Grid % per_x) < abs(lx) ) then
        lx = lx - Grid % per_x
      else if( abs(lx + Grid % per_x) < abs(lx) ) then
        lx = lx + Grid % per_x
      end if
      if( abs(ly - Grid % per_y) < abs(ly) ) then
        ly = ly - Grid % per_y
      else if( abs(ly + Grid % per_y) < abs(ly) ) then
        ly = ly + Grid % per_y
      end if
      if( abs(lz - Grid % per_z) < abs(lz) ) then
        lz = lz - Grid % per_z
      else if( abs(lz + Grid % per_z) < abs(lz) ) then
        lz = lz + Grid % per_z
      end if

      Grid % weight_c2n(i_cell, n) = 1.0 + lambda_x * lx   &
                                         + lambda_y * ly   &
                                         + lambda_z * lz

    end do  ! through cells surrounding the node

  end do  ! through nodes

  !---------------------------!
  !   Normalize the weights   !
  !---------------------------!
  do n = 1, Grid % n_nodes

    ! Add total weights for all cells
    tot = 0.0
    do i_cell = 1, Grid % nodes_n_cells(n)
      tot = tot + Grid % weight_c2n(i_cell, n)
    end do

    ! Divide each weight with total
    do i_cell = 1, Grid % nodes_n_cells(n)
      Grid % weight_c2n(i_cell, n) = Grid % weight_c2n(i_cell, n) / tot
    end do

  end do

  end subroutine

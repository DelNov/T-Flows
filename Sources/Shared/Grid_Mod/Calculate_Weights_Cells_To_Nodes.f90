!==============================================================================!
  subroutine Grid_Mod_Calculate_Weights_Cells_To_Nodes(grid)
!------------------------------------------------------------------------------!
!   Computes weights for interpolation from cell to nodes                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, n, i_cell
  logical :: here
  real    :: lx, ly, lz, rx, ry, rz, lambda_x, lambda_y, lambda_z
  real    :: ixx, iyy, izz, ixz, iyz, ixy, d
  real    :: a11, a12, a13, a21, a22, a23, a31, a32, a33
  real    :: tot, weights_sorted(64)
!==============================================================================!

  ! Allocate memory
  allocate(grid % nodes_weight_c(size(grid % nodes_c,1),size(grid % nodes_c,2)))
  grid % nodes_weight_c = 0.0

  !-------------------------------------------!
  !   Browse through all nodes to calculate   !
  !     weights from cells surrounding it     !
  !-------------------------------------------!
  do n = 1, grid % n_nodes
    rx = 0.0; ry = 0.0; rz = 0.0
    ixx = 0.0; iyy = 0.0; izz = 0.0; ixz = 0.0; iyz = 0.0; ixy = 0.0

    ! Loop on cells surrounding the node
    do i_cell = 1, grid % nodes_n_cells(n)
      c = grid % nodes_c(i_cell, n)

      ! Default distances
      lx = grid % xc(c) - grid % xn(n)
      ly = grid % yc(c) - grid % yn(n)
      lz = grid % zc(c) - grid % zn(n)

      ! Check if periodicity in all three directions and correct distances
      if( abs(lx - grid % per_x) < abs(lx) ) then
        lx = lx - grid % per_x
      else if( abs(lx + grid % per_x) < abs(lx) ) then
        lx = lx + grid % per_x
      end if
      if( abs(ly - grid % per_y) < abs(ly) ) then
        ly = ly - grid % per_y
      else if( abs(ly + grid % per_y) < abs(ly) ) then
        ly = ly + grid % per_y
      end if
      if( abs(lz - grid % per_z) < abs(lz) ) then
        lz = lz - grid % per_z
      else if( abs(lz + grid % per_z) < abs(lz) ) then
        lz = lz + grid % per_z
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

    do i_cell = 1, grid % nodes_n_cells(n)
      c  = grid % nodes_c(i_cell, n)

      ! Default distances
      lx = grid % xc(c) - grid % xn(n)
      ly = grid % yc(c) - grid % yn(n)
      lz = grid % zc(c) - grid % zn(n)

      ! Check if periodicity in all three directions and correct distances
      if( abs(lx - grid % per_x) < abs(lx) ) then
        lx = lx - grid % per_x
      else if( abs(lx + grid % per_x) < abs(lx) ) then
        lx = lx + grid % per_x
      end if
      if( abs(ly - grid % per_y) < abs(ly) ) then
        ly = ly - grid % per_y
      else if( abs(ly + grid % per_y) < abs(ly) ) then
        ly = ly + grid % per_y
      end if
      if( abs(lz - grid % per_z) < abs(lz) ) then
        lz = lz - grid % per_z
      else if( abs(lz + grid % per_z) < abs(lz) ) then
        lz = lz + grid % per_z
      end if

      grid % nodes_weight_c(i_cell, n) = 1.0             &
                                       + lambda_x * lx   &
                                       + lambda_y * ly   &
                                       + lambda_z * lz

    end do  ! through cells surrounding the node

  end do  ! through nodes

  !---------------------------!
  !   Normalize the weights   !
  !---------------------------!
  do n = 1, grid % n_nodes

    ! Add total weights for all cells
    tot = 0.0
    do i_cell = 1, grid % nodes_n_cells(n)
      tot = tot + grid % nodes_weight_c(i_cell, n)
    end do

    ! Divide each weight with total
    do i_cell = 1, grid % nodes_n_cells(n)
      grid % nodes_weight_c(i_cell, n) = grid % nodes_weight_c(i_cell, n) / tot
    end do

  end do

! ! Debugging
! write(200 + this_proc, '(a)')  'List of nodes with their cells weights'
! do n = 1, grid % n_nodes
!   here = .false.
!   do i_cell = 1, grid % nodes_n_cells(n)
!     c = grid % nodes_c(i_cell, n)
!     if(grid % comm % cell_proc(c) .eq. this_proc) here = .true.
!   end do
!   if(here) then
!     weights_sorted(1:grid % nodes_n_cells(n)) =  &
!     grid % nodes_weight_c(1:grid % nodes_n_cells(n), n)
!     call Sort_Mod_Real(weights_sorted(1:grid % nodes_n_cells(n)))
!     write(200 + this_proc, '(i7.7, i3, 99f7.4)')  &
!                 grid % comm % node_glo(n),        &
!                 grid % nodes_n_cells(n),          &
!                 weights_sorted(1:grid % nodes_n_cells(n))
!   end if
! end do

  end subroutine

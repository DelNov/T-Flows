!==============================================================================!
  subroutine Field_Mod_Calculate_Grad_Matrix_Cells_To_Nodes(flow)
!------------------------------------------------------------------------------!
!   Calculates gradient matrix.                                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: c, i_cel, n
  real                     :: dx, dy, dz
  real                     :: jac, g_inv(6)
!==============================================================================!

  ! Take alias
  grid => flow % pnt_grid

  ! Innitialize all entries to zero
  do n = 1, grid % n_nodes
    flow % grad_c2n(1,n) = 0.0
    flow % grad_c2n(2,n) = 0.0
    flow % grad_c2n(3,n) = 0.0
    flow % grad_c2n(4,n) = 0.0
    flow % grad_c2n(5,n) = 0.0
    flow % grad_c2n(6,n) = 0.0
  end do

  !---------------------------------!
  !   Compute the gradient matrix   !
  !---------------------------------!
  do n = 1, grid % n_nodes

    ! Browse through node's cells
    do i_cel = 1, grid % nodes_n_cells(n)
      c  = grid % nodes_c(i_cel, n)
      dx = grid % xc(c) - grid % xn(n)
      dy = grid % yc(c) - grid % yn(n)
      dz = grid % zc(c) - grid % zn(n)

      flow % grad_c2n(1,n) = flow % grad_c2n(1,n) + dx * dx    ! 1,1
      flow % grad_c2n(2,n) = flow % grad_c2n(2,n) + dy * dy    ! 2,2
      flow % grad_c2n(3,n) = flow % grad_c2n(3,n) + dz * dz    ! 3,3
      flow % grad_c2n(4,n) = flow % grad_c2n(4,n) + dx * dy    ! 1,2  &  2,1
      flow % grad_c2n(5,n) = flow % grad_c2n(5,n) + dx * dz    ! 1,3  &  3,1
      flow % grad_c2n(6,n) = flow % grad_c2n(6,n) + dy * dz    ! 2,3  &  3,2
    end do

  end do

  !--------------------------------!
  !   Find the inverse of matrix   !
  !--------------------------------!
  do n = 1, grid % n_nodes
    jac = flow % grad_c2n(1,n) * flow % grad_c2n(2,n) * flow % grad_c2n(3,n)  &
        - flow % grad_c2n(1,n) * flow % grad_c2n(6,n) * flow % grad_c2n(6,n)  &
        - flow % grad_c2n(4,n) * flow % grad_c2n(4,n) * flow % grad_c2n(3,n)  &
        + flow % grad_c2n(4,n) * flow % grad_c2n(5,n) * flow % grad_c2n(6,n)  &
        + flow % grad_c2n(4,n) * flow % grad_c2n(5,n) * flow % grad_c2n(6,n)  &
        - flow % grad_c2n(5,n) * flow % grad_c2n(5,n) * flow % grad_c2n(2,n)

    g_inv(1) = +(  flow % grad_c2n(2,n) * flow % grad_c2n(3,n)  &
                 - flow % grad_c2n(6,n) * flow % grad_c2n(6,n) ) / (jac+TINY)
    g_inv(2) = +(  flow % grad_c2n(1,n) * flow % grad_c2n(3,n)  &
                 - flow % grad_c2n(5,n) * flow % grad_c2n(5,n) ) / (jac+TINY)
    g_inv(3) = +(  flow % grad_c2n(1,n) * flow % grad_c2n(2,n)  &
                 - flow % grad_c2n(4,n) * flow % grad_c2n(4,n) ) / (jac+TINY)
    g_inv(4) = -(  flow % grad_c2n(4,n) * flow % grad_c2n(3,n)  &
                 - flow % grad_c2n(5,n) * flow % grad_c2n(6,n) ) / (jac+TINY)
    g_inv(5) = +(  flow % grad_c2n(4,n) * flow % grad_c2n(6,n)  &
                 - flow % grad_c2n(5,n) * flow % grad_c2n(2,n) ) / (jac+TINY)
    g_inv(6) = -(  flow % grad_c2n(1,n) * flow % grad_c2n(6,n)  &
                 - flow % grad_c2n(4,n) * flow % grad_c2n(5,n) ) / (jac+TINY)

    flow % grad_c2n(1,n) = g_inv(1)
    flow % grad_c2n(2,n) = g_inv(2)
    flow % grad_c2n(3,n) = g_inv(3)
    flow % grad_c2n(4,n) = g_inv(4)
    flow % grad_c2n(5,n) = g_inv(5)
    flow % grad_c2n(6,n) = g_inv(6)
  end do

  end subroutine

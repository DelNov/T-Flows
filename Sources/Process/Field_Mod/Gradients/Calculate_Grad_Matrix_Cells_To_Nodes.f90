!==============================================================================!
  subroutine Calculate_Grad_Matrix_Cells_To_Nodes(Flow)
!------------------------------------------------------------------------------!
!   Calculates gradient matrix.                                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type) :: Flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: c, i_cel, n
  real                     :: dx, dy, dz
  real                     :: jac, g_inv(6)
!==============================================================================!

  ! Take alias
  grid => Flow % pnt_grid

  !--------------------------------------------!
  !   Initialize gradient matrices for nodes   !
  !--------------------------------------------!
  do n = 1, grid % n_nodes
    Flow % grad_c2n(1,n) = 0.0
    Flow % grad_c2n(2,n) = 0.0
    Flow % grad_c2n(3,n) = 0.0
    Flow % grad_c2n(4,n) = 0.0
    Flow % grad_c2n(5,n) = 0.0
    Flow % grad_c2n(6,n) = 0.0
  end do

  !-----------------------------------------------!
  !   Compute the gradient matrix for all nodes   !
  !-----------------------------------------------!
  do n = 1, grid % n_nodes

    ! Browse through node's cells
    do i_cel = 1, grid % nodes_n_cells(n)
      c  = grid % nodes_c(i_cel, n)
      dx = grid % xc(c) - grid % xn(n)
      dy = grid % yc(c) - grid % yn(n)
      dz = grid % zc(c) - grid % zn(n)

      Flow % grad_c2n(1,n) = Flow % grad_c2n(1,n) + dx * dx    ! 1,1
      Flow % grad_c2n(2,n) = Flow % grad_c2n(2,n) + dy * dy    ! 2,2
      Flow % grad_c2n(3,n) = Flow % grad_c2n(3,n) + dz * dz    ! 3,3
      Flow % grad_c2n(4,n) = Flow % grad_c2n(4,n) + dx * dy    ! 1,2  &  2,1
      Flow % grad_c2n(5,n) = Flow % grad_c2n(5,n) + dx * dz    ! 1,3  &  3,1
      Flow % grad_c2n(6,n) = Flow % grad_c2n(6,n) + dy * dz    ! 2,3  &  3,2
    end do

  end do

  !--------------------------------!
  !   Find the inverse of matrix   !
  !--------------------------------!
  do n = 1, grid % n_nodes
    jac = Flow % grad_c2n(1,n) * Flow % grad_c2n(2,n) * Flow % grad_c2n(3,n)  &
        - Flow % grad_c2n(1,n) * Flow % grad_c2n(6,n) * Flow % grad_c2n(6,n)  &
        - Flow % grad_c2n(4,n) * Flow % grad_c2n(4,n) * Flow % grad_c2n(3,n)  &
        + Flow % grad_c2n(4,n) * Flow % grad_c2n(5,n) * Flow % grad_c2n(6,n)  &
        + Flow % grad_c2n(4,n) * Flow % grad_c2n(5,n) * Flow % grad_c2n(6,n)  &
        - Flow % grad_c2n(5,n) * Flow % grad_c2n(5,n) * Flow % grad_c2n(2,n)

    g_inv(1) = +(  Flow % grad_c2n(2,n) * Flow % grad_c2n(3,n)  &
                 - Flow % grad_c2n(6,n) * Flow % grad_c2n(6,n) ) / (jac+TINY)
    g_inv(2) = +(  Flow % grad_c2n(1,n) * Flow % grad_c2n(3,n)  &
                 - Flow % grad_c2n(5,n) * Flow % grad_c2n(5,n) ) / (jac+TINY)
    g_inv(3) = +(  Flow % grad_c2n(1,n) * Flow % grad_c2n(2,n)  &
                 - Flow % grad_c2n(4,n) * Flow % grad_c2n(4,n) ) / (jac+TINY)
    g_inv(4) = -(  Flow % grad_c2n(4,n) * Flow % grad_c2n(3,n)  &
                 - Flow % grad_c2n(5,n) * Flow % grad_c2n(6,n) ) / (jac+TINY)
    g_inv(5) = +(  Flow % grad_c2n(4,n) * Flow % grad_c2n(6,n)  &
                 - Flow % grad_c2n(5,n) * Flow % grad_c2n(2,n) ) / (jac+TINY)
    g_inv(6) = -(  Flow % grad_c2n(1,n) * Flow % grad_c2n(6,n)  &
                 - Flow % grad_c2n(4,n) * Flow % grad_c2n(5,n) ) / (jac+TINY)

    Flow % grad_c2n(1,n) = g_inv(1)
    Flow % grad_c2n(2,n) = g_inv(2)
    Flow % grad_c2n(3,n) = g_inv(3)
    Flow % grad_c2n(4,n) = g_inv(4)
    Flow % grad_c2n(5,n) = g_inv(5)
    Flow % grad_c2n(6,n) = g_inv(6)
  end do

  end subroutine

!==============================================================================!
  subroutine Calculate_Grad_Matrix_Nodes_To_Cells(Flow)
!------------------------------------------------------------------------------!
!   Calculates gradient matrix.                                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type) :: Flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: c, i_nod, n
  real                     :: dx, dy, dz
  real                     :: jac, g_inv(6)
!==============================================================================!

  ! Take alias
  grid => Flow % pnt_grid

  !--------------------------------------------!
  !   Initialize gradient matrices for cells   !
  !--------------------------------------------!
  do c = 1, grid % n_cells
    Flow % grad_n2c(1,c) = 0.0
    Flow % grad_n2c(2,c) = 0.0
    Flow % grad_n2c(3,c) = 0.0
    Flow % grad_n2c(4,c) = 0.0
    Flow % grad_n2c(5,c) = 0.0
    Flow % grad_n2c(6,c) = 0.0
  end do

  !-----------------------------------------------!
  !   Compute the gradient matrix for all nodes   !
  !-----------------------------------------------!
  do c = 1, grid % n_cells

    ! Browse through cell's nodes
    do i_nod = 1, abs(grid % cells_n_nodes(c))
      n  = grid % cells_n(i_nod, c)
      dx = grid % xn(n) - grid % xc(c)
      dy = grid % yn(n) - grid % yc(c)
      dz = grid % zn(n) - grid % zc(c)

      Flow % grad_n2c(1,c) = Flow % grad_n2c(1,c) + dx * dx    ! 1,1
      Flow % grad_n2c(2,c) = Flow % grad_n2c(2,c) + dy * dy    ! 2,2
      Flow % grad_n2c(3,c) = Flow % grad_n2c(3,c) + dz * dz    ! 3,3
      Flow % grad_n2c(4,c) = Flow % grad_n2c(4,c) + dx * dy    ! 1,2  &  2,1
      Flow % grad_n2c(5,c) = Flow % grad_n2c(5,c) + dx * dz    ! 1,3  &  3,1
      Flow % grad_n2c(6,c) = Flow % grad_n2c(6,c) + dy * dz    ! 2,3  &  3,2
    end do

  end do

  !--------------------------------!
  !   Find the inverse of matrix   !
  !--------------------------------!
  do c = 1, grid % n_cells
    jac = Flow % grad_n2c(1,c) * Flow % grad_n2c(2,c) * Flow % grad_n2c(3,c)  &
        - Flow % grad_n2c(1,c) * Flow % grad_n2c(6,c) * Flow % grad_n2c(6,c)  &
        - Flow % grad_n2c(4,c) * Flow % grad_n2c(4,c) * Flow % grad_n2c(3,c)  &
        + Flow % grad_n2c(4,c) * Flow % grad_n2c(5,c) * Flow % grad_n2c(6,c)  &
        + Flow % grad_n2c(4,c) * Flow % grad_n2c(5,c) * Flow % grad_n2c(6,c)  &
        - Flow % grad_n2c(5,c) * Flow % grad_n2c(5,c) * Flow % grad_n2c(2,c)

    g_inv(1) = +(  Flow % grad_n2c(2,c) * Flow % grad_n2c(3,c)  &
                 - Flow % grad_n2c(6,c) * Flow % grad_n2c(6,c) ) / (jac+TINY)
    g_inv(2) = +(  Flow % grad_n2c(1,c) * Flow % grad_n2c(3,c)  &
                 - Flow % grad_n2c(5,c) * Flow % grad_n2c(5,c) ) / (jac+TINY)
    g_inv(3) = +(  Flow % grad_n2c(1,c) * Flow % grad_n2c(2,c)  &
                 - Flow % grad_n2c(4,c) * Flow % grad_n2c(4,c) ) / (jac+TINY)
    g_inv(4) = -(  Flow % grad_n2c(4,c) * Flow % grad_n2c(3,c)  &
                 - Flow % grad_n2c(5,c) * Flow % grad_n2c(6,c) ) / (jac+TINY)
    g_inv(5) = +(  Flow % grad_n2c(4,c) * Flow % grad_n2c(6,c)  &
                 - Flow % grad_n2c(5,c) * Flow % grad_n2c(2,c) ) / (jac+TINY)
    g_inv(6) = -(  Flow % grad_n2c(1,c) * Flow % grad_n2c(6,c)  &
                 - Flow % grad_n2c(4,c) * Flow % grad_n2c(5,c) ) / (jac+TINY)

    Flow % grad_n2c(1,c) = g_inv(1)
    Flow % grad_n2c(2,c) = g_inv(2)
    Flow % grad_n2c(3,c) = g_inv(3)
    Flow % grad_n2c(4,c) = g_inv(4)
    Flow % grad_n2c(5,c) = g_inv(5)
    Flow % grad_n2c(6,c) = g_inv(6)
  end do

  end subroutine

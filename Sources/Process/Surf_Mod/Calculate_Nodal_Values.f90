!==============================================================================!
  subroutine Surf_Mod_Calculate_Nodal_Values(surf, phi)
!------------------------------------------------------------------------------!
!   Calculates nodal values of variable phi                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: phi_n => r_node_01,  &  ! value at the (static) grid nodes
                      ncn   => i_node_01      ! number of cells around nodes
!------------------------------------------------------------------------------!
!   Be careful with the above variables from Work_Mod.  They are used by       !
!   two subroutines in Surf_Mod, hence values shouldn't be changed elsewhere.  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
  type(Var_Type),  target :: phi
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: c, ni, n
  real                     :: xc, yc, zc
!==============================================================================!

  ! Take aliases
  grid => surf % pnt_grid

  !-------------------------------------------------------!
  !   Find number of inside cells surrounding each node   !
  !-------------------------------------------------------!
  do n = 1, grid % n_nodes
    ncn(n) = 0
  end do
  do c = 1, grid % n_cells
    do ni = 1, grid % cells_n_nodes(c)  ! local node number
      n = grid % cells_n(ni, c)         ! global node number

      ! Increase number of cells surrounding the this node by one
      ncn(n) = ncn(n) + 1
    end do
  end do

  !-------------------------------------------------!
  !   Compute values of variable phi on the nodes   !
  !-------------------------------------------------!

  do n = 1, grid % n_nodes
    phi_n(n) = 0.0
  end do

  ! Add them up for all nodes from cells ...
  do c = 1, grid % n_cells ! think about it: - grid % comm % n_buff_cells
    xc = grid % xc(c)
    yc = grid % yc(c)
    zc = grid % zc(c)
    do ni = 1, grid % cells_n_nodes(c)
      n = grid % cells_n(ni, c)
      phi_n(n) = phi_n(n) + phi % n(c)
    end do
  end do

  ! ... and divide by number of cells surrounding each node
  do n = 1, grid % n_nodes
    phi_n(n) = phi_n(n) / ncn(n)
  end do

  end subroutine

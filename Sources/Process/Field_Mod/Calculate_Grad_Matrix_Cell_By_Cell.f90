!==============================================================================!
  subroutine Field_Mod_Calculate_Grad_Matrix_Cell_By_Cell(flow)
!------------------------------------------------------------------------------!
!   Calculates gradient matrix for cells, browsing through cells.              !
!                                                                              !
!   Essentially a front-end for Field_Mod_Calculate_Grad_Matrix_For_Cell       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: c, s, i_fac
  real                     :: dx_c, dy_c, dz_c
  real                     :: jac, g_inv(6)
!==============================================================================!

  ! Take alias
  grid => flow % pnt_grid

  !-----------------------------------------------!
  !   Calculate gradient matrices for all cells   !
  !-----------------------------------------------!
  do c = 1, grid % n_cells
    call Field_Mod_Calculate_Grad_Matrix_For_Cell(flow, c)
  end do

  end subroutine

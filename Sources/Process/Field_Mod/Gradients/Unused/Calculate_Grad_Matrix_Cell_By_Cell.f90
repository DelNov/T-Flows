!==============================================================================!
  subroutine Calculate_Grad_Matrix_Cell_By_Cell(Flow)
!------------------------------------------------------------------------------!
!   Calculates gradient matrix for cells, browsing through cells.              !
!                                                                              !
!   Essentially a front-end for Field_Mod_Calculate_Grad_Matrix_For_Cell       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type) :: Flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: c, s, i_fac
  real                     :: dx_c, dy_c, dz_c
  real                     :: jac, g_inv(6)
!==============================================================================!

  ! Take alias
  Grid => Flow % pnt_grid

  !-----------------------------------------------!
  !   Calculate gradient matrices for all cells   !
  !-----------------------------------------------!
  do c = 1, Grid % n_cells
    call Flow % Calculate_Grad_Matrix_For_Cell(c)
  end do

  end subroutine

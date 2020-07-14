!==============================================================================!
  subroutine Multiphase_Mod_Vof_Interpolate_Faces_Cells(grid,                 &
                                                        var_face, var_cell)
!------------------------------------------------------------------------------!
!   Interpolates a variable from nodes to cell center                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: var_face(1:grid % n_faces)
  real            :: var_cell(-grid % n_bnd_cells:grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer :: n, s, c, c1, c2, i_fac
  real    :: sum1, sum2
!==============================================================================!

  do c = 1, grid % n_cells
    sum1 = 0.0; sum2 = 0.0
    do i_fac = 1, grid % cells_n_faces(c)
      s = grid % cells_f(i_fac, c)
      sum1 = sum1 + grid % weight_faces(i_fac, c) * var_face(s)
      sum2 = sum2 + grid % weight_faces(i_fac, c)
    end do
    var_cell(c) = sum1 / sum2
  end do

  end subroutine

!==============================================================================!
  subroutine Initialize_New_Numbers(Grid)
!------------------------------------------------------------------------------!
!>  The subroutine is designed to set up reference mappings for nodes, cells,
!>  and faces within a Grid object. It ensures that each element points to its
!>  own index, useful when saving entire grids (not sub-domains) in
!>  pre-processing stages.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid  !! grid under consideration
!-----------------------------------[Locals]-----------------------------------!
  integer :: n, c, s
!==============================================================================!

  do n = 1, Grid % n_nodes
    Grid % new_n(n) = n
  end do

  do c = -Grid % n_bnd_cells, Grid % n_cells
    Grid % new_c(c) = c
    Grid % old_c(c) = c
  end do

  do s = 1, Grid % n_faces + Grid % n_shadows
    Grid % new_f(s) = s
    Grid % old_f(s) = s
  end do

  end subroutine

!==============================================================================!
  subroutine Initialize_New_Numbers(Grid)
!------------------------------------------------------------------------------!
!   Set (almost) all new_ and old_ arrays to point to their own selves         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
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

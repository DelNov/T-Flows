!==============================================================================!
  subroutine Grid_Mod_Initialize_New_Numbers(grid)
!------------------------------------------------------------------------------!
!   Set (almost) all new_ and old_ arrays to point to their own selves         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: n, c, s
!==============================================================================!

  do n = 1, grid % n_nodes
    grid % new_n(n) = n
  end do

  do c = -grid % n_bnd_cells, grid % n_cells
    grid % new_c(c) = c
    grid % old_c(c) = c
  end do

  do s = 1, grid % n_faces + grid % n_shadows
    grid % new_f(s) = s
    grid % old_f(s) = s
  end do

  end subroutine

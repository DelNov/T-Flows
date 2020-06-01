!==============================================================================!
  subroutine Grid_Mod_Allocate_New_Numbers(grid, nn, nb, nc, nf)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: nn, nb, nc, nf
!==============================================================================!

  allocate (grid % new_n(    nn));  grid % new_n(:) = 0
  allocate (grid % new_c(-nb:nc));  grid % new_c(:) = 0
  allocate (grid % new_f(    nf));  grid % new_f(:) = 0

  allocate (grid % old_c(-nb:nc));  grid % old_c(:) = 0

  end subroutine

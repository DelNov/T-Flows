!==============================================================================!
  integer function N_Cells_In_Bnd_Color(grid, bc)
!------------------------------------------------------------------------------!
!   Counts and marks (with new_c) boundary cells in the given boundary color   !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: bc
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, cnt
!==============================================================================!

  ! Nullify on entry
  cnt = 0
  grid % new_c(:) = 0

  do c = -grid % n_bnd_cells, -1
    if( grid % bnd_cond % color(c) .eq. bc ) then
      cnt = cnt + 1
      grid % new_c(c) = cnt
    end if
  end do

  ! Return a sum of all marked boundary cells
  N_Cells_In_Bnd_Color = cnt

  end function


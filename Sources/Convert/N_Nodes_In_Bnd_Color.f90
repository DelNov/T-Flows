!==============================================================================!
  integer function N_Nodes_In_Bnd_Color(grid, bc)
!------------------------------------------------------------------------------!
!   Counts and marks (with new_n) nodes in the given boundary color            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: bc
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, i_nod, n, cnt
!==============================================================================!

  ! Nullify on entry
  cnt = 0
  grid % new_n(:) = 0

  ! Browse through all boundary cells and mark nodes
  ! of those cells in the given boundary condition
  do c = -grid % n_bnd_cells, -1
    if( grid % bnd_cond % color(c) .eq. bc ) then
      do i_nod = 1, grid % cells_n_nodes(c)
        n = grid % cells_n(i_nod, c)
        if(grid % new_n(n) .eq. 0) then  ! hasn't been marked yet
          cnt = cnt + 1
          grid % new_n(n) = cnt
        end if
      end do
    end if
  end do

  ! Return a sum of all marked nodes
  N_Nodes_In_Bnd_Color = cnt

  end function


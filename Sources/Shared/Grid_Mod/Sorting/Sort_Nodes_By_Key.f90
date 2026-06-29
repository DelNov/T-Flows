!==============================================================================!
  subroutine Sort_Nodes_By_Key(Grid, key)
!------------------------------------------------------------------------------!
!>  Sort nodes by the intiger array passed with argument "key".
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type), intent(inout) :: Grid                 !! parent grid
  integer                         :: key(Grid % n_nodes)  !! index array
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, c, n, i_nod
!==============================================================================!

  Assert(PROGRAM_NAME .eq. 'Convert')

  Assert(allocated(Grid % xn))

  write(*,'(a)', advance='no') ' # Sorting node coordinates ...'

  !----------------------------!
  !   Store old node numbers   !
  !----------------------------!
  do n = 1, Grid % n_nodes
    Grid % old_n(n) = n
  end do

  !--------------------------------------------------!
  !   Sort new numbers according to three criteria   !
  !--------------------------------------------------!
  call Sort % Int_Carry_Int(key, Grid % old_n(1:Grid % n_nodes))

  ! This is a bit of a bluff
  do n = 1, Grid % n_nodes
    Grid % new_n(Grid % old_n(n)) = n
  end do

  ! Sort node coordinates
  call Sort % Real_By_Index(Grid % n_nodes, Grid % xn, Grid % new_n)
  call Sort % Real_By_Index(Grid % n_nodes, Grid % yn, Grid % new_n)
  call Sort % Real_By_Index(Grid % n_nodes, Grid % zn, Grid % new_n)

  !-----------------------------------------------!
  !   Do the sorting of data pertinent to cells   !
  !-----------------------------------------------!
  do c = -Grid % n_bnd_cells, Grid % n_cells
    do i_nod = 1, abs(Grid % cells_n_nodes(c))
      n = Grid % cells_n(i_nod, c)
      Grid % cells_n(i_nod, c) = Grid % new_n(n)
    end do
  end do

  !-----------------------------------------------!
  !   Do the sorting of data pertinent to faces   !
  !-----------------------------------------------!
  do s = 1, Grid % n_faces + Grid % n_shadows
    do i_nod = 1, Grid % faces_n_nodes(s)
      n = Grid % faces_n(i_nod, s)
      Grid % faces_n(i_nod, s) = Grid % new_n(n)
    end do
  end do

  print '(a)', ' done!'

  end subroutine

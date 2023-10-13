!==============================================================================!
  subroutine Sort_Nodes_By_Coordinates(Grid)
!------------------------------------------------------------------------------!
!   Sorts nodes by their geometrical positions.  Used by Insert_Buildings      !
!   from Convert, before placing a city on a ground (STL) file.                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type), intent(inout) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, c, n, i_nod
!==============================================================================!

  Assert(PROGRAM_NAME .eq. 'Convert')

  Assert(allocated(Grid % xn))

  write(*,'(a)', advance='no') ' # Sorting node coordinates ...'

  ! It is higlhy unlikely this error will ever occur
  if(Math % Approx_Real(maxval(Grid % xn(:)), 0.0) .and.  &
     Math % Approx_Real(maxval(Grid % yn(:)), 0.0) .and.  &
     Math % Approx_Real(maxval(Grid % zn(:)), 0.0) .and.  &
     Math % Approx_Real(minval(Grid % xn(:)), 0.0) .and.  &
     Math % Approx_Real(minval(Grid % yn(:)), 0.0) .and.  &
     Math % Approx_Real(minval(Grid % zn(:)), 0.0)) then
    call Message % Error(64,                                                   &
              'You called the function Grid % Sort_Nodes_By_Coordinates  ' //  &
              'before the node coordinates have been read.  Something is ' //  &
              'wrong with logic of the algorithm, you better re-asses it.',    &
              file=__FILE__, line=__LINE__)
  end if

  !----------------------------!
  !   Store old node numbers   !
  !----------------------------!
  do n = 1, Grid % n_nodes
    Grid % old_n(n) = n
  end do

  !--------------------------------------------------!
  !   Sort new numbers according to three criteria   !
  !--------------------------------------------------!
  call Sort % Three_Real_Carry_Int(Grid % xn(1:Grid % n_nodes),  &
                                   Grid % yn(1:Grid % n_nodes),  &
                                   Grid % zn(1:Grid % n_nodes),  &
                                   Grid % old_n(1:Grid % n_nodes))
  ! This is a bit of a bluff
  do n = 1, Grid % n_nodes
    Grid % new_n(Grid % old_n(n)) = n
  end do

  !-----------------------------------------------!
  !   Do the sorting of data pertinent to cells   !
  !-----------------------------------------------!
  do c = -Grid % n_bnd_cells, Grid % n_cells
    do i_nod = 1, Grid % cells_n_nodes(c)
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

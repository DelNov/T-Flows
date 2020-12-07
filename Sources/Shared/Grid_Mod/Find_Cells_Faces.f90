!==============================================================================!
  subroutine Grid_Mod_Find_Cells_Faces(grid)
!------------------------------------------------------------------------------!
!   Find faces surrounding each cell; including boundary cells                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s   ! counters
!==============================================================================!

  if(any(grid % cells_n_faces(:) .ne. 0)) then
    print *, '# WARNING: Seems you are looking for cells'' faces' //  &
             ' although this information has already been found!'
  end if

  !---------------------------------------!
  !   Browse through all faces to form:   !
  !        cells_n_faces, cells_f         !
  !---------------------------------------!
  grid % cells_n_faces(:) = 0
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(2, s)

    grid % cells_n_faces(c1) = grid % cells_n_faces(c1) + 1
    grid % cells_f(grid % cells_n_faces(c1), c1) = s

    grid % cells_n_faces(c2) = grid % cells_n_faces(c2) + 1
    grid % cells_f(grid % cells_n_faces(c2), c2) = s
  end do

  end subroutine

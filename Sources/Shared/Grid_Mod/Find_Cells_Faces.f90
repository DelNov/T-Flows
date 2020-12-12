!==============================================================================!
  subroutine Grid_Mod_Find_Cells_Faces(grid)
!------------------------------------------------------------------------------!
!   Find faces surrounding each cell; including boundary cells                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, i_cel, s   ! counters
!==============================================================================!

  if(any(grid % cells_n_faces(:) .ne. 0)) then
    print *, '# NOTE: Seems you are looking for cells'' faces'
    print *, '# although this information has already been found!'
    print *, '# No harm done, just a little note.'
  end if

  grid % cells_n_faces(:) = 0

  do s = 1, grid % n_faces
    do i_cel = 1, 2
      c = grid % faces_c(i_cel, s)

      grid % cells_n_faces(c) = grid % cells_n_faces(c) + 1
      grid % cells_f(grid % cells_n_faces(c), c) = s
    end do
  end do

  end subroutine

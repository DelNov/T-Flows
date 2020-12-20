!==============================================================================!
  subroutine Grid_Mod_Find_Cells_Faces(grid)
!------------------------------------------------------------------------------!
!   Find faces surrounding each cell; including boundary cells                 !
!                                                                              !
!   Please note that it takes shadow faces into account as well                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, i_cel, s, sh   ! counters
  real    :: dist_s, dist_sh
!==============================================================================!

  if(any(grid % cells_n_faces(:) .ne. 0)) then
    print *, '# NOTE: Seems you are looking for cells'' faces'
    print *, '# although this information has already been found!'
    print *, '# No harm done, just a little note.'
  end if

  grid % cells_n_faces(:) = 0

  do s = 1, grid % n_faces
    do i_cel = 1, 2
      c = grid % faces_c(i_cel, s)  ! would be c1 and c2 in most of the code

      grid % cells_n_faces(c) = grid % cells_n_faces(c) + 1
      grid % cells_f(grid % cells_n_faces(c), c) = s

      sh = grid % faces_s(s)        ! take the shadow face
      if(sh > 0) then
        dist_s  = Math_Mod_Distance(                                &
                       grid % xc(c),  grid % yc(c),  grid % zc(c),  &
                       grid % xf(s),  grid % yf(s),  grid % zf(s))
        dist_sh = Math_Mod_Distance(                                &
                       grid % xc(c),  grid % yc(c),  grid % zc(c),  &
                       grid % xf(sh), grid % yf(sh), grid % zf(sh))
        if(dist_sh < dist_s) then
          grid % cells_f(grid % cells_n_faces(c), c) = s
        end if
      end if
    end do
  end do

  end subroutine

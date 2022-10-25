!==============================================================================!
  subroutine Find_Cells_Faces(Grid)
!------------------------------------------------------------------------------!
!   Find faces surrounding each cell; including boundary cells                 !
!                                                                              !
!   Please note that it takes shadow faces into account as well                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, i_cel, s, sh   ! counters
  real    :: dist_s, dist_sh
!==============================================================================!

  if(any(Grid % cells_n_faces(:) .ne. 0)) then
    print *, '# NOTE: Seems you are looking for cells'' faces'
    print *, '# although this information has already been found!'
    print *, '# No harm done, just a little note.'
  end if

  Grid % cells_n_faces(:) = 0

  do s = 1, Grid % n_faces
    do i_cel = 1, 2
      c = Grid % faces_c(i_cel, s)  ! would be c1 and c2 in most of the code

      Grid % cells_n_faces(c) = Grid % cells_n_faces(c) + 1
      call Adjust_First_Dim(Grid % cells_n_faces(c), Grid % cells_f)
      Grid % cells_f(Grid % cells_n_faces(c), c) = s

      sh = Grid % faces_s(s)        ! take the shadow face
      if(sh > 0) then
        dist_s  = Math % Distance(                                  &
                       Grid % xc(c),  Grid % yc(c),  Grid % zc(c),  &
                       Grid % xf(s),  Grid % yf(s),  Grid % zf(s))
        dist_sh = Math % Distance(                                  &
                       Grid % xc(c),  Grid % yc(c),  Grid % zc(c),  &
                       Grid % xf(sh), Grid % yf(sh), Grid % zf(sh))
        if(dist_sh < dist_s) then
          Grid % cells_f(Grid % cells_n_faces(c), c) = s
        end if
      end if
    end do
  end do

  end subroutine

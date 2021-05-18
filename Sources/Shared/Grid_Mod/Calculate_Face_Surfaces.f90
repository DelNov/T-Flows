!==============================================================================!
  subroutine Calculate_Face_Surfaces(Grid)
!------------------------------------------------------------------------------!
!   Calculate the face surface areas                                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, i_nod, j_nod, n
  real    :: xt(MAX_FACES_N_NODES), yt(MAX_FACES_N_NODES), zt(MAX_FACES_N_NODES)
!==============================================================================!

  do s = 1, Grid % n_faces

    ! Copy face node coordinates to a local array for easier handling
    do i_nod = 1, Grid % faces_n_nodes(s)  ! local node counter
      n = Grid % faces_n(i_nod, s)         ! global node number
      xt(i_nod) = Grid % xn(n)
      yt(i_nod) = Grid % yn(n)
      zt(i_nod) = Grid % zn(n)
    end do

    ! Cell face components
    Grid % sx(s) = 0.0
    Grid % sy(s) = 0.0
    Grid % sz(s) = 0.0
    do i_nod = 1, Grid % faces_n_nodes(s)
      j_nod = i_nod + 1
      if(j_nod > Grid % faces_n_nodes(s)) j_nod = 1
      Grid % sx(s) = Grid % sx(s)  &
                   + (yt(j_nod) - yt(i_nod)) * (zt(j_nod) + zt(i_nod))
      Grid % sy(s) = Grid % sy(s)  &
                   + (zt(j_nod) - zt(i_nod)) * (xt(j_nod) + xt(i_nod))
      Grid % sz(s) = Grid % sz(s)  &
                   + (xt(j_nod) - xt(i_nod)) * (yt(j_nod) + yt(i_nod))
    end do
    Grid % sx(s) = 0.5 * Grid % sx(s)
    Grid % sy(s) = 0.5 * Grid % sy(s)
    Grid % sz(s) = 0.5 * Grid % sz(s)

  end do

  print *, '# Face surfaces calculated !'

  end subroutine

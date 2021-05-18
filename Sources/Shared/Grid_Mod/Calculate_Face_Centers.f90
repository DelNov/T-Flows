!==============================================================================!
  subroutine Calculate_Face_Centers(Grid)
!------------------------------------------------------------------------------!
!   Calculate the face centers from nodal coordinates.                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, i_nod, n, nn
  real    :: xt(MAX_FACES_N_NODES), yt(MAX_FACES_N_NODES), zt(MAX_FACES_N_NODES)
!==============================================================================!

  do s = 1, Grid % n_faces

    ! Copy face node coordinates to a local array for easier handling
    do i_nod = 1, Grid % faces_n_nodes(s)  ! local node counter
      n = Grid % faces_n(i_nod ,s)  ! global node number
      xt(i_nod) = Grid % xn(n)
      yt(i_nod) = Grid % yn(n)
      zt(i_nod) = Grid % zn(n)
    end do

    ! Barycenters
    nn = Grid % faces_n_nodes(s)
    Grid % xf(s) = sum( xt(1:nn) ) / real(nn)
    Grid % yf(s) = sum( yt(1:nn) ) / real(nn)
    Grid % zf(s) = sum( zt(1:nn) ) / real(nn)

  end do

  print *, '# Face centers calculated !'

  end subroutine

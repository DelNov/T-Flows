!==============================================================================!
  subroutine Grid_Mod_Calculate_Face_Centers(grid)
!------------------------------------------------------------------------------!
!   Calculate the face centers from nodal coordinates.                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, i_nod, n, nn
  real    :: xt(MAX_FACES_N_NODES), yt(MAX_FACES_N_NODES), zt(MAX_FACES_N_NODES)
!==============================================================================!

  do s = 1, grid % n_faces

    ! Copy face node coordinates to a local array for easier handling
    do i_nod = 1, grid % faces_n_nodes(s)  ! local node counter
      n = grid % faces_n(i_nod ,s)  ! global node number
      xt(i_nod) = grid % xn(n)
      yt(i_nod) = grid % yn(n)
      zt(i_nod) = grid % zn(n)
    end do

    ! Barycenters
    nn = grid % faces_n_nodes(s)
    grid % xf(s) = sum( xt(1:nn) ) / real(nn)
    grid % yf(s) = sum( yt(1:nn) ) / real(nn)
    grid % zf(s) = sum( zt(1:nn) ) / real(nn)

  end do

  print *, '# Face centers calculated !'

  end subroutine

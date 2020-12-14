!==============================================================================!
  subroutine Grid_Mod_Calculate_Face_Surfaces(grid)
!------------------------------------------------------------------------------!
!   Calculate the face surface areas                                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, i_nod, j_nod, n
  real    :: xt(MAX_FACES_N_NODES), yt(MAX_FACES_N_NODES), zt(MAX_FACES_N_NODES)
!==============================================================================!

  do s = 1, grid % n_faces

    ! Copy face node coordinates to a local array for easier handling
    do i_nod = 1, grid % faces_n_nodes(s)  ! local node counter
      n = grid % faces_n(i_nod, s)         ! global node number
      xt(i_nod) = grid % xn(n)
      yt(i_nod) = grid % yn(n)
      zt(i_nod) = grid % zn(n)
    end do

    ! Cell face components
    grid % sx(s) = 0.0
    grid % sy(s) = 0.0
    grid % sz(s) = 0.0
    do i_nod = 1, grid % faces_n_nodes(s)
      j_nod = i_nod + 1
      if(j_nod > grid % faces_n_nodes(s)) j_nod = 1
      grid % sx(s) = grid % sx(s)  &
                   + (yt(j_nod) - yt(i_nod)) * (zt(j_nod) + zt(i_nod))
      grid % sy(s) = grid % sy(s)  &
                   + (zt(j_nod) - zt(i_nod)) * (xt(j_nod) + xt(i_nod))
      grid % sz(s) = grid % sz(s)  &
                   + (xt(j_nod) - xt(i_nod)) * (yt(j_nod) + yt(i_nod))
    end do
    grid % sx(s) = 0.5 * grid % sx(s)
    grid % sy(s) = 0.5 * grid % sy(s)
    grid % sz(s) = 0.5 * grid % sz(s)

  end do

  print *, '# Face surfaces calculated !'

  end subroutine

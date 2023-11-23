!==============================================================================!
  subroutine Calculate_Face_Centers(Grid)
!------------------------------------------------------------------------------!
!>  Calculate the face centers from nodal coordinates.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid  !! grid under consideration
!-----------------------------------[Locals]-----------------------------------!
  integer           :: s, i_nod, m, n, nn
  real, allocatable :: xf(:), yf(:), zf(:)
!==============================================================================!

  ! Allocate memory for face's node coordinates
  m = size(Grid % faces_n, 1)
  allocate(xf(m))
  allocate(yf(m))
  allocate(zf(m))

  ! Do the actual calculation
  do s = 1, Grid % n_faces

    ! Copy face node coordinates to a local array for easier handling
    do i_nod = 1, Grid % faces_n_nodes(s)  ! local node counter
      n = Grid % faces_n(i_nod ,s)  ! global node number
      xf(i_nod) = Grid % xn(n)
      yf(i_nod) = Grid % yn(n)
      zf(i_nod) = Grid % zn(n)
    end do

    ! Barycenters
    nn = Grid % faces_n_nodes(s)
    Grid % xf(s) = sum( xf(1:nn) ) / real(nn)
    Grid % yf(s) = sum( yf(1:nn) ) / real(nn)
    Grid % zf(s) = sum( zf(1:nn) ) / real(nn)

  end do

  print *, '# Face centers calculated !'

  end subroutine

!==============================================================================!
  subroutine Faces_Surface(Grid, s, sx, sy, sz)
!------------------------------------------------------------------------------!
!>  Calculate the face surface area for one face.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)     :: Grid        !! grid under consideration
  integer, intent(in)  :: s           !! face number
  real,    intent(out) :: sx, sy, sz  !! component of the face surface area
!-----------------------------------[Locals]-----------------------------------!
  integer           :: i_nod, j_nod, m, n
  real, allocatable :: xf(:), yf(:), zf(:)
!==============================================================================!

  ! Allocate memory for face's node coordinates
  m = Grid % faces_n_nodes(s)
  allocate(xf(m))
  allocate(yf(m))
  allocate(zf(m))

  ! Copy face node coordinates to a local array for easier handling
  do i_nod = 1, Grid % faces_n_nodes(s)  ! local node counter
    n = Grid % faces_n(i_nod, s)         ! global node number
    xf(i_nod) = Grid % xn(n)
    yf(i_nod) = Grid % yn(n)
    zf(i_nod) = Grid % zn(n)
  end do

  ! Cell face components
  sx = 0.0
  sy = 0.0
  sz = 0.0
  do i_nod = 1, Grid % faces_n_nodes(s)
    j_nod = i_nod + 1
    if(j_nod > Grid % faces_n_nodes(s)) j_nod = 1
    sx = sx + (yf(j_nod) - yf(i_nod)) * (zf(j_nod) + zf(i_nod))
    sy = sy + (zf(j_nod) - zf(i_nod)) * (xf(j_nod) + xf(i_nod))
    sz = sz + (xf(j_nod) - xf(i_nod)) * (yf(j_nod) + yf(i_nod))
  end do
  sx = 0.5 * sx
  sy = 0.5 * sy
  sz = 0.5 * sz

  end subroutine

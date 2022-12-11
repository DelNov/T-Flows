!==============================================================================!
  subroutine Calculate_Face_Surfaces(Grid)
!------------------------------------------------------------------------------!
!   Calculate the face surface areas                                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c1, c2, s, i_nod, j_nod, m, n
  real              :: dx, dy, dz
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
      n = Grid % faces_n(i_nod, s)         ! global node number
      xf(i_nod) = Grid % xn(n)
      yf(i_nod) = Grid % yn(n)
      zf(i_nod) = Grid % zn(n)
    end do

    ! Cell face components
    Grid % sx(s) = 0.0
    Grid % sy(s) = 0.0
    Grid % sz(s) = 0.0

    ! I am afraid that this algorithm computes surfaces with wrong sign, but
    ! it might work afterall since the nodes in "fn" structures (HEX, TET,
    ! WED and PYR in Convert_Mod) also have nodes sorted in the direction
    ! which points inwards a cell and not outwards
    do i_nod = 1, Grid % faces_n_nodes(s)
      j_nod = i_nod + 1
      if(j_nod > Grid % faces_n_nodes(s)) j_nod = 1
      Grid % sx(s) = Grid % sx(s)  &
                   + (yf(j_nod) - yf(i_nod)) * (zf(j_nod) + zf(i_nod))
      Grid % sy(s) = Grid % sy(s)  &
                   + (zf(j_nod) - zf(i_nod)) * (xf(j_nod) + xf(i_nod))
      Grid % sz(s) = Grid % sz(s)  &
                   + (xf(j_nod) - xf(i_nod)) * (yf(j_nod) + yf(i_nod))
    end do
    Grid % sx(s) = 0.5 * Grid % sx(s)
    Grid % sy(s) = 0.5 * Grid % sy(s)
    Grid % sz(s) = 0.5 * Grid % sz(s)

    ! Perform assertion for inside faces
    ! (Boundary cells and faces are not
    !  formed at this point yet anyhow.)
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)
    if(c2 .gt. 0) then
      dx = grid % xc(c2) - grid % xc(c1)
      dy = grid % yc(c2) - grid % yc(c1)
      dz = grid % zc(c2) - grid % zc(c1)
      Assert(Grid % sx(s) * dx + Grid % sy(s) * dy + Grid % sz(s) * dz > 0.0)
    end if

  end do

  print *, '# Face surfaces calculated !'

  end subroutine

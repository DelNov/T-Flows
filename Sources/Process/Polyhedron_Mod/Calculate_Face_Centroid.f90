!==============================================================================!
  subroutine Calculate_Face_Centroid(Polyhedron, face, xf, yf, zf)
!------------------------------------------------------------------------------!
!   Calculate face's centroid, this information is not stored in Polyhedron    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Polyhedron
  integer, intent(in)    :: face
  real,    intent(out)   :: xf, yf, zf
!-----------------------------------[Locals]-----------------------------------!
  integer :: i_nod, i
!==============================================================================!

  xf = 0.0
  yf = 0.0
  zf = 0.0
  do i_nod = 1, Polyhedron % faces_n_nodes(face)
    i  = Polyhedron % faces_n(face, i_nod)
    xf = xf + Polyhedron % nodes_xyz(i, 1)
    yf = yf + Polyhedron % nodes_xyz(i, 2)
    zf = zf + Polyhedron % nodes_xyz(i, 3)
  end do
  xf = xf / real(Polyhedron % faces_n_nodes(face))
  yf = yf / real(Polyhedron % faces_n_nodes(face))
  zf = zf / real(Polyhedron % faces_n_nodes(face))

  end subroutine

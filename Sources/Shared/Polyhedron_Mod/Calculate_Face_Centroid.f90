!==============================================================================!
  subroutine Calculate_Face_Centroid(Pol, face, xf, yf, zf)
!------------------------------------------------------------------------------!
!   Calculate face's centroid, this information is not stored in Polyhedron    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Pol
  integer, intent(in)    :: face
  real,    intent(out)   :: xf, yf, zf
!-----------------------------------[Locals]-----------------------------------!
  integer :: i_nod, i
!==============================================================================!

  xf = 0.0
  yf = 0.0
  zf = 0.0
  do i_nod = 1, Pol % faces_n_nodes(face)
    i  = Pol % faces_n(face, i_nod)
    xf = xf + Pol % nodes_xyz(i, 1)
    yf = yf + Pol % nodes_xyz(i, 2)
    zf = zf + Pol % nodes_xyz(i, 3)
  end do
  if(Pol % faces_n_nodes(face) > 0) then
    xf = xf / real(Pol % faces_n_nodes(face))
    yf = yf / real(Pol % faces_n_nodes(face))
    zf = zf / real(Pol % faces_n_nodes(face))
  end if

  end subroutine

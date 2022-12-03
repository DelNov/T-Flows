!==============================================================================!
  function Facets_Vert_Coords(Stl, fac, v)
!------------------------------------------------------------------------------!
!   Returns cooridnates of facets's vertex                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Stl_Type)     :: Stl
  real, dimension(3)  :: Facets_Vert_Coords
  integer, intent(in) :: fac
  integer, intent(in) :: v
!==============================================================================!

  Assert(v <= 3)

  Facets_Vert_Coords(1) = Stl % xn(Stl % cells_n(v, -fac))
  Facets_Vert_Coords(2) = Stl % yn(Stl % cells_n(v, -fac))
  Facets_Vert_Coords(3) = Stl % zn(Stl % cells_n(v, -fac))

  end function

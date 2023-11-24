!==============================================================================!
  function Facets_Vert_Coords(Stl, fac, v)
!------------------------------------------------------------------------------!
!>  This function, belonging to the Stl_Mod module, returns the coordinates
!>  of a specified vertex (numbered locally) in a specified facet.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Stl_Type)     :: Stl                 !! parent Stl_Type object
  real, dimension(3)  :: Facets_Vert_Coords  !! output, vertex coordinates
  integer, intent(in) :: fac                 !! facet number
  integer, intent(in) :: v                   !! facet's local vertex
!==============================================================================!

  Assert(v <= 3)

  Facets_Vert_Coords(1) = Stl % xn(Stl % cells_n(v, -fac))
  Facets_Vert_Coords(2) = Stl % yn(Stl % cells_n(v, -fac))
  Facets_Vert_Coords(3) = Stl % zn(Stl % cells_n(v, -fac))

  end function

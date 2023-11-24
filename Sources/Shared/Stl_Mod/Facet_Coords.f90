!==============================================================================!
  function Facet_Coords(Stl, fac)
!------------------------------------------------------------------------------!
!>  This function, belonging to the Stl_Mod module, calculates and returns the
!>  coordinates of the center of gravity (centroid) of a specified facet in an
!>  Stl object.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Stl_Type)     :: Stl           !! parent Stl_Type object
  real, dimension(3)  :: Facet_Coords  !! output, centroid of the facet
  integer, intent(in) :: fac           !! facet number
!==============================================================================!

  Facet_Coords(1) = (  Stl % xn(Stl % cells_n(1,-fac))  &
                     + Stl % xn(Stl % cells_n(2,-fac))  &
                     + Stl % xn(Stl % cells_n(3,-fac))) * ONE_THIRD

  Facet_Coords(2) = (  Stl % yn(Stl % cells_n(1,-fac))  &
                     + Stl % yn(Stl % cells_n(2,-fac))  &
                     + Stl % yn(Stl % cells_n(3,-fac))) * ONE_THIRD

  Facet_Coords(3) = (  Stl % zn(Stl % cells_n(1,-fac))  &
                     + Stl % zn(Stl % cells_n(2,-fac))  &
                     + Stl % zn(Stl % cells_n(3,-fac))) * ONE_THIRD

  end function

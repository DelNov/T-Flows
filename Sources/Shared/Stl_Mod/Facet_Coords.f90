!==============================================================================!
  function Facet_Coords(Stl, fac)
!------------------------------------------------------------------------------!
!   Returns facet cooridnates, facet's center of gravity, in other words.      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Stl_Type)     :: Stl
  real, dimension(3)  :: Facet_Coords
  integer, intent(in) :: fac
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

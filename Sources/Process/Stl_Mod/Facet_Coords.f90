!==============================================================================!
  function Facet_Coords(Stl, fac)
!------------------------------------------------------------------------------!
!   Returns facet cooridnates                                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Stl_Type)     :: Stl
  real, dimension(3)  :: Facet_Coords
  integer, intent(in) :: fac
!==============================================================================!

  Facet_Coords(1) = (Stl % x(1,fac) + Stl % x(2,fac) + Stl % x(3,fac)) / 3.0
  Facet_Coords(2) = (Stl % y(1,fac) + Stl % y(2,fac) + Stl % y(3,fac)) / 3.0
  Facet_Coords(3) = (Stl % z(1,fac) + Stl % z(2,fac) + Stl % z(3,fac)) / 3.0

  end function

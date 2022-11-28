!==============================================================================!
  function Facet_Normal(Stl, fac)
!------------------------------------------------------------------------------!
!   Returns facet cooridnates                                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Stl_Type)     :: Stl
  real, dimension(3)  :: Facet_Normal
  integer, intent(in) :: fac
!==============================================================================!

  Facet_Normal(1) = Stl % nx(fac)
  Facet_Normal(2) = Stl % ny(fac)
  Facet_Normal(3) = Stl % nz(fac)

  end function

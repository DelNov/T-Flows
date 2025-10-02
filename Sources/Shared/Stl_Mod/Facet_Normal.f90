!==============================================================================!
  function Facet_Normal(Stl, fac)
!------------------------------------------------------------------------------!
!>  This function, belonging to the Stl_Mod module, returns the facet's normal.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Stl_Type)     :: Stl           !! parent Stl_Type object
  real, dimension(3)  :: Facet_Normal  !! output, facet's normal vector
  integer, intent(in) :: fac           !! facet number
!==============================================================================!

  Facet_Normal(1) = Stl % nx(-fac)
  Facet_Normal(2) = Stl % ny(-fac)
  Facet_Normal(3) = Stl % nz(-fac)

  end function

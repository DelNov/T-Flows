!==============================================================================!
  integer function Facet_Body(Stl, fac)
!------------------------------------------------------------------------------!
!>  This function, part of the Stl_Mod module, returns the body identifier
!>  associated with a given facet in an STL object. It is a straightforward
!>  lookup function that retrieves the body number for a specific facet.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Stl_Type)     :: Stl  !! parent Stl_Type object
  integer, intent(in) :: fac  !! facet number
!==============================================================================!

  Facet_Body = Stl % body_c(-fac)

  end function

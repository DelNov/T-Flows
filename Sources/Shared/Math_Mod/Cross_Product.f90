!==============================================================================!
  function Cross_Product(Math, a, b)
!------------------------------------------------------------------------------!
!   This is a prototype of a small module which would contain some basic       !
!   mathematic and related functions.                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type)               :: Math
  real, dimension(3)             :: Cross_Product
  real, dimension(3), intent(in) :: a, b
!==============================================================================!

  Cross_Product(1) = a(2) * b(3) - a(3) * b(2)
  Cross_Product(2) = a(3) * b(1) - a(1) * b(3)
  Cross_Product(3) = a(1) * b(2) - a(2) * b(1)

  end function

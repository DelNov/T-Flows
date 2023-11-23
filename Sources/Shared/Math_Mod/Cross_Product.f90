!==============================================================================!
  pure function Cross_Product(Math, a, b)
!------------------------------------------------------------------------------!
!>  Calculates a cross (vectorial) product of two three-dimensional vectors.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type),   intent(in) :: Math           !! parent class
  real, dimension(3)             :: Cross_Product  !! resulting vector
  real, dimension(3), intent(in) :: a, b           !! operand of the product
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Math)
!==============================================================================!

  Cross_Product(1) = a(2) * b(3) - a(3) * b(2)
  Cross_Product(2) = a(3) * b(1) - a(1) * b(3)
  Cross_Product(3) = a(1) * b(2) - a(2) * b(1)

  end function

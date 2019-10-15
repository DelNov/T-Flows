!==============================================================================!
  function Math_Mod_Cross_Product(a, b)
!------------------------------------------------------------------------------!
!   This is a prototype of a small module which would contain some basic       !
!   mathematic and related functions.                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, dimension(3)             :: Math_Mod_Cross_Product
  real, dimension(3), intent(in) :: a, b
!==============================================================================!

  Math_Mod_Cross_Product(1) = a(2) * b(3) - a(3) * b(2)
  Math_Mod_Cross_Product(2) = a(3) * b(1) - a(1) * b(3)
  Math_Mod_Cross_Product(3) = a(1) * b(2) - a(2) * b(1)

  end function

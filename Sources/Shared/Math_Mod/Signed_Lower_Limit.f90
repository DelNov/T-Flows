!==============================================================================!
  real function Signed_Lower_Limit(Math, a, limit)
!------------------------------------------------------------------------------!
!   Calculates a lower limit, preserving a sign                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type) :: Math
  real             :: a
  real, intent(in) :: limit  ! lower limit, usually a small number
!==============================================================================!

  ! Assume it won't change
  Signed_Lower_Limit = a

  ! If smaller than the limit
  if(abs(a) < limit) then
    Signed_Lower_Limit = sign(limit, a)
  end if

  end function

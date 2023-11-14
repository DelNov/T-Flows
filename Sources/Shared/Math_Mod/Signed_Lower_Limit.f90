!==============================================================================!
  pure real function Signed_Lower_Limit(Math, a, limit)
!------------------------------------------------------------------------------!
!>  Calculates a lower limit (of a small number), preserving a sign.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type), intent(in) :: Math   !! parent class
  real,             intent(in) :: a      !! input value, usually a small number
  real,             intent(in) :: limit  !! lower limit
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Math)
!==============================================================================!

  ! Assume it won't change
  Signed_Lower_Limit = a

  ! If smaller than the limit
  if(abs(a) < limit) then
    Signed_Lower_Limit = sign(limit, a)
  end if

  end function

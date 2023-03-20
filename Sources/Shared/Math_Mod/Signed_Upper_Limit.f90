!==============================================================================!
  pure real function Signed_Upper_Limit(Math, a, limit)
!------------------------------------------------------------------------------!
!   Calculates an upper limit, preserving the sign                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type), intent(in) :: Math
  real,             intent(in) :: a
  real,             intent(in) :: limit  ! upper limit, usually a big number
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Math)
!==============================================================================!

  ! Assume it won't change
  Signed_Upper_Limit = a

  ! If smaller than the limit
  if(abs(a) > limit) then
    Signed_Upper_Limit = sign(limit, a)
  end if

  end function

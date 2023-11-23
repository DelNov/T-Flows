!==============================================================================!
  pure logical function Smaller_Real(Math, a, b, tol)
!------------------------------------------------------------------------------!
!>  Returns .true. if a is approximatelly smaller than b, .false. otherwise.
!------------------------------------------------------------------------------!
!   The approximation is controlled by optional parameter "tol".  If it is     !
!   not given, the function will use DEFAULT_TOLERANCE specified in Math_Mod   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type), intent(in) :: Math  !! parent class
  real,             intent(in) :: a, b  !! operands
  real, optional,   intent(in) :: tol   !! tolerance
!-----------------------------------[Locals]-----------------------------------!
  real :: tolerance
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Math)
!==============================================================================!

  if( .not. present(tol) ) then
    tolerance = DEFAULT_TOLERANCE
  else
    tolerance = tol
  end if

  if( a < b - tolerance ) then
    Smaller_Real = .true.
  else
    Smaller_Real = .false.
  end if

  end function

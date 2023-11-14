!==============================================================================!
  pure logical function Approx_Real(Math, a, b, tol)
!------------------------------------------------------------------------------!
!>  Returns .true. if a is approximatelly equal to "b", .false. otherwise.
!>  The approximation is controlled by optional parameter "tol".  If it is
!>  not given, the function will use DEFAULT_TOLERANCE specified in Math_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type), intent(in) :: Math  !! parent class
  real,             intent(in) :: a, b  !! numbers being compared
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

  if( abs(a - b) <= tolerance ) then
    Approx_Real = .true.
  else
    Approx_Real = .false.
  end if

  end function

!==============================================================================!
  logical function Math_Mod_Approx_Real(a, b, tol)
!------------------------------------------------------------------------------!
!   Returns .true. if a is approximatelly equal to "b", .false. otherwise.     !
!                                                                              !
!   The approximation is controlled by optional parameter "tol".  If it is     !
!   not given, the function will use DEFAULT_TOLERANCE specified in Math_Mod   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real           :: a, b
  real, optional :: tol
!-----------------------------------[Locals]-----------------------------------!
  real :: tolerance
!==============================================================================!

  if( .not. present(tol) ) then
    tolerance = DEFAULT_TOLERANCE
  else
    tolerance = tol
  end if

  if( abs(a - b) <= tolerance ) then
    Math_Mod_Approx_Real = .true.
  else
    Math_Mod_Approx_Real = .false.
  end if

  end function

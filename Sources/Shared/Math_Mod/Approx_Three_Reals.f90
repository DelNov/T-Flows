!==============================================================================!
  logical function Approx_Three_Reals(Math, a, b, tol)
!------------------------------------------------------------------------------!
!   Returns .true. if a(1:3) is approximatelly equal to b(1:3), and .false.    !
!   if it is not the case                                                      !
!                                                                              !
!   The approximation is controlled by optional parameter "tol".  If it is     !
!   not given, the function will use DEFAULT_TOLERANCE specified in Math_Mod   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type) :: Math
  real             :: a(1:3), b(1:3)
  real, optional   :: tol
!-----------------------------------[Locals]-----------------------------------!
  real :: tolerance
!==============================================================================!

  if( .not. present(tol) ) then
    tolerance = DEFAULT_TOLERANCE
  else
    tolerance = tol
  end if

  Approx_Three_Reals = Math % Approx_Real(a(1), b(1), tolerance) .and.  &
                       Math % Approx_Real(a(2), b(2), tolerance) .and.  &
                       Math % Approx_Real(a(3), b(3), tolerance)

  end function

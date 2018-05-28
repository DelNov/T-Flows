!==============================================================================!
  logical function Approx(a, b, tol)
!------------------------------------------------------------------------------!
!   Returns .true. if a is approximatelly equal to "b", .false. otherwise.     !
!                                                                              !
!   The approximation is controlled by optional parameter "tol".  If it is     !
!   not provided, the function will use value 1.0e-9 as a default.             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real           :: a, b
  real, optional :: tol
!-----------------------------------[Locals]-----------------------------------!
  real :: tolerance
!==============================================================================!

  if( .not. present(tol) ) then
    tolerance = 1.e-9
  else
    tolerance = tol
  end if

  if( (a < (b + tolerance)) .and. (a > (b - tolerance)) ) then
    approx = .true.
  else
    approx = .false.
  end if

  end function Approx

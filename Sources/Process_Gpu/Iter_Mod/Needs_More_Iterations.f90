!==============================================================================!
  logical function Needs_More_Iterations(Iter, Flow, n_dom)
!------------------------------------------------------------------------------!
!>  Determines if additional iterations are required and returns .true. if
!>  that is the case, or .false. if no further iterations are needed.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Iter_Type), intent(inout) :: Iter      !! parent, singleton object Iter
  type(Field_Type), intent(in)    :: Flow(MD)  !! all flow domains in simulation
  integer,          intent(in)    :: n_dom     !! actual number of domains
!==============================================================================!

  !--------------------------------!
  !   Increase current iteration   !
  !--------------------------------!
  Iter % current_iteration = Iter % current_iteration + 1

  !------------------------------------------------------------------!
  !   Minimum number of iterations hasn't been reached -> continue   !
  !------------------------------------------------------------------!
  if(Iter % Current() .le. Iter % Get_Min()) then
    Needs_More_Iterations = .true.
    return
  end if

  !-----------------------------------------------------------!
  !   Maximum number of iterations has been reached -> stop   !
  !-----------------------------------------------------------!
  if(Iter % Current() .gt. Iter % Get_Max()) then
    Iter % current_iteration = 0
    Needs_More_Iterations = .false.
    return
  end if

  !--------------------------------------------!
  !   You are in between minimum and maximum   !
  !   number of iterations, check residuals    !
  !--------------------------------------------!

  ! Residuals reached the tolerance level -> stop
  if(Iter % Max_Fields_Residual(Flow, n_dom) <= Iter % Get_Tol()) then
    Iter % current_iteration = 0
    Needs_More_Iterations = .false.
    return

  ! Residuals did not reached the tolerance level -> continue
  else
    Needs_More_Iterations = .true.
    return
  end if

  ! You shouldn't be here
  Assert(.false.)

  end function

!==============================================================================!
  subroutine Iterations(Rc)
!------------------------------------------------------------------------------!
!   Reads control file for variables which control outer iteration loop        !
!---------------------------------[Arguments]----------------------------------!
  class(Read_Controls_Type), intent(in) :: Rc
!-----------------------------------[Locals]-----------------------------------!
  integer :: max_out     ! max number of inner iterations
  integer :: min_out     ! min number of inner iterations
  real    :: simple_tol  ! tolerance for SIMPLE algorithm
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Rc)
!==============================================================================!

  call Control % Max_Simple_Iterations(max_out)
  call Control % Min_Simple_Iterations(min_out)
  call Control % Tolerance_For_Simple_Algorithm(simple_tol)

  call Iter % Set_Max(max_out)
  call Iter % Set_Min(min_out)
  call Iter % Set_Tol(simple_tol)

  end subroutine

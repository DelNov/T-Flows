!==============================================================================!
  subroutine Iterations(Rc)
!------------------------------------------------------------------------------!
!>  This is s a simplified version from the same subroutine in Process_Cpu
!>  as it reads only boundary conditions releated to momentum and enthalpy
!>  conservation equations.  Hopefully, as more modules are ported to
!>  Process_Gpu, this source file will get closer and closer to its origin
!>  from Process_Cpu.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Read_Controls_Type), intent(in) :: Rc  !! parent class
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


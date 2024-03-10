!==============================================================================!
  subroutine Native_Solvers(Rc, Flow)
!------------------------------------------------------------------------------!
!>  This is s a simplified version from the same subroutine in Process_Cpu
!>  as it reads only boundary conditions releated to momentum and enthalpy
!>  conservation equations.  Hopefully, as more modules are ported to
!>  Process_Gpu, this source file will get closer and closer to its origin
!>  from Process_Cpu.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Read_Controls_Type), intent(in) :: Rc    !! parent class
  type(Field_Type), target              :: Flow  !! flow object
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: ui, phi
  integer                  :: i, sc
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Rc)
!==============================================================================!

  ! Give some sign
  if(First_Proc())  &
    print '(a)', ' # Reading info about native (home-brewed) solvers'

  ! Take alias
  Grid => Flow % pnt_grid

  !-------------------------!
  !   Related to momentum   !
  !-------------------------!
  do i = 1, 3
    if(i .eq. 1) ui => Flow % u
    if(i .eq. 2) ui => Flow % v
    if(i .eq. 3) ui => Flow % w
    call Control % Tolerance_For_Momentum_Solver     (ui % tol)
    call Control % Max_Iterations_For_Momentum_Solver(ui % miter)
  end do

  !-------------------------!
  !   Related to pressure   !
  !-------------------------!
  call Control % Tolerance_For_Pressure_Solver     (Flow % pp % tol)
  call Control % Max_Iterations_For_Pressure_Solver(Flow % pp % miter)

  !------------------------------!
  !   Related to heat transfer   !
  !------------------------------!
  if(Flow % heat_transfer) then
    call Control % Tolerance_For_Energy_Solver      (Flow % t % tol)
    call Control % Max_Iterations_For_Energy_Solver (Flow % t % miter)
  end if

  !--------------------------------!
  !   Related to passive scalars   !
  !--------------------------------!
  if(Flow % n_scalars .gt. 0) then

    ! Although there is a slight assymetry here (allocations are done in
    ! Create_Field), this call is needed here to reserve memory to store
    ! info on linear solvers. Fields are still allocated in Create_Field
    allocate(Flow % scalar(Flow % n_scalars))

    ! Now, with basic allocation done, it is safe to read solver options
    do sc = 1, Flow % n_scalars
      phi => Flow % scalar(sc)
      call Control % Tolerance_For_Scalars_Solver     (phi % tol)
      call Control % Max_Iterations_For_Scalars_Solver(phi % miter)
    end do

  end if

  end subroutine

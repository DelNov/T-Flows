!==============================================================================!
  subroutine Native_Solvers(Rc, Flow, Turb, Vof)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to read and set up the parameters for various
!>  native (home-brewed) solvers from the control file.  It configures solvers
!>  for different aspects of the flow simulation such as momentum, pressure,
!>  wall distance, potential, heat transfer, multiphase flow, passive scalars,
!>  and turbulence quantities.
!------------------------------------------------------------------------------!
!   Note on good practice                                                      !
!                                                                              !
!   * Default values, outlined in Documents/all_control_keywords, should be    !
!     defined only in Control_Mod, not be scattered around the code.  In other !
!     words, Control_Mod changes less frequently than other parts of the code, !
!     so it is safer to place default values there.                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Read_Controls_Type), intent(in) :: Rc    !! parent class
  type(Field_Type), target              :: Flow  !! flow object
  type(Turb_Type),  target              :: Turb  !! turbulence object
  type(Vof_Type),   target              :: Vof   !! VOF object
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: tq, ui, phi
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
    call Control % Solver_For_Momentum               (ui % solver)
    call Control % Preconditioner_For_System_Matrix  (ui % prec)
    call Control % Tolerance_For_Momentum_Solver     (ui % tol)
    call Control % Max_Iterations_For_Momentum_Solver(ui % miter)
  end do

  !-------------------------!
  !   Related to pressure   !
  !-------------------------!
  call Control % Solver_For_Pressure               (Flow % pp % solver)
  call Control % Preconditioner_For_System_Matrix  (Flow % pp % prec)
  call Control % Tolerance_For_Pressure_Solver     (Flow % pp % tol)
  call Control % Max_Iterations_For_Pressure_Solver(Flow % pp % miter)

  !------------------------------!
  !   Related to wall distance   !
  !------------------------------!
  call Control % Solver_For_Wall_Distance            (Flow % wall_dist % solver)
  call Control % Preconditioner_For_System_Matrix    (Flow % wall_dist % prec)
  call Control % Tolerance_For_Wall_Distance_Solver  (Flow % wall_dist % tol)
  call Control % Max_Iterations_For_Wall_Distance_Solver  &
                                                     (Flow % wall_dist % miter)

  !--------------------------!
  !   Related to potential   !  (for flow field initialization)
  !--------------------------!
  call Control % Solver_For_Potential               (Flow % pot % solver)
  call Control % Preconditioner_For_System_Matrix   (Flow % pot % prec)
  call Control % Tolerance_For_Potential_Solver     (Flow % pot % tol)
  call Control % Max_Iterations_For_Potential_Solver(Flow % pot % miter)

  !------------------------------!
  !   Related to heat transfer   !
  !------------------------------!
  if(Flow % heat_transfer) then
    call Control % Solver_For_Energy                (Flow % t % solver)
    call Control % Preconditioner_For_System_Matrix (Flow % t % prec)
    call Control % Tolerance_For_Energy_Solver      (Flow % t % tol)
    call Control % Max_Iterations_For_Energy_Solver (Flow % t % miter)
  end if

  !--------------------------------!
  !   Related to multiphase Flow   !
  !--------------------------------!
  if(Flow % with_interface) then
    call Control % Solver_For_Vof                  (Vof % fun % solver)
    call Control % Preconditioner_For_System_Matrix(Vof % fun % prec)
    call Control % Tolerance_For_Vof_Solver        (Vof % fun % tol)
    call Control % Max_Iterations_For_Vof_Solver   (Vof % fun % miter)
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
      call Control % Solver_For_Scalars               (phi % solver)
      call Control % Preconditioner_For_System_Matrix (phi % prec)
      call Control % Tolerance_For_Scalars_Solver     (phi % tol)
      call Control % Max_Iterations_For_Scalars_Solver(phi % miter)
    end do

  end if

  !------------------------------!
  !   All turbuelnt quantities   !
  !------------------------------!
  do i = 1, 12
    if(i .eq.  1) tq => Turb % kin
    if(i .eq.  2) tq => Turb % eps
    if(i .eq.  3) tq => Turb % zeta
    if(i .eq.  4) tq => Turb % f22
    if(i .eq.  5) tq => Turb % vis
    if(i .eq.  6) tq => Turb % t2
    if(i .eq.  7) tq => Turb % uu
    if(i .eq.  8) tq => Turb % vv
    if(i .eq.  9) tq => Turb % ww
    if(i .eq. 10) tq => Turb % uv
    if(i .eq. 11) tq => Turb % uw
    if(i .eq. 12) tq => Turb % vw
    call Control % Solver_For_Turbulence               (tq % solver)
    call Control % Preconditioner_For_System_Matrix    (tq % prec)
    call Control % Tolerance_For_Turbulence_Solver     (tq % tol)
    call Control % Max_Iterations_For_Turbulence_Solver(tq % miter)
  end do

  end subroutine

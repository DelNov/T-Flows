!==============================================================================!
  subroutine Native_Solvers(Rc, Flow, turb, Vof, Sol)
!------------------------------------------------------------------------------!
!   Reads details about native solvers from control file.                      !
!                                                                              !
!   Good practice: default values, outlined in Documents/all_control_keywords, !
!   should be defined only in Control_Mod, not be scattered around the code.   !
!   In other words, Control_Mod changes less frequently than other parts of    !
!   the code, so it is safer to place default values there.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Read_Controls_Type)  :: Rc
  type(Field_Type),   target :: Flow
  type(Turb_Type),    target :: turb
  type(Vof_Type),     target :: Vof
  type(Solver_Type),  target :: Sol
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: tq, ui, phi
  character(SL)            :: name
  integer                  :: i, sc
!==============================================================================!

  ! Take alias
  Grid => Flow % pnt_grid

  !----------------------------------!
  !   Gradient computation methods   !
  !----------------------------------!

  ! Tolerance and max iterations for computation of gradients with Gauss method
  call Control_Mod_Tolerance_For_Gauss_Gradients(Flow % gauss_tol,    .false.)
  call Control_Mod_Max_Gauss_Gradients_Iterations(Flow % gauss_miter, .false.)
  call Control_Mod_Max_Least_Squares_Gradients_Iterations(Flow % least_miter,  &
                                                          .false.)

  !-------------------------!
  !   Related to momentum   !
  !-------------------------!
  do i = 1, 3
    if(i .eq. 1) ui => Flow % u
    if(i .eq. 2) ui => Flow % v
    if(i .eq. 3) ui => Flow % w
    call Control_Mod_Solver_For_Momentum               (ui % solver)
    call Control_Mod_Preconditioner_For_System_Matrix  (ui % prec)
    call Control_Mod_Tolerance_For_Momentum_Solver     (ui % tol)
    call Control_Mod_Max_Iterations_For_Momentum_Solver(ui % mniter)
  end do

  !-------------------------!
  !   Related to pressure   !
  !-------------------------!
  call Control_Mod_Solver_For_Pressure               (Flow % pp % solver)
  call Control_Mod_Preconditioner_For_System_Matrix  (Flow % pp % prec)
  call Control_Mod_Tolerance_For_Pressure_Solver     (Flow % pp % tol)
  call Control_Mod_Max_Iterations_For_Pressure_Solver(Flow % pp % mniter)

  !------------------------------!
  !   Related to wall distance   !
  !------------------------------!
  call Control_Mod_Solver_For_Wall_Distance                 &
                                                     (Flow % wall_dist % solver)
  call Control_Mod_Preconditioner_For_System_Matrix         &
                                                     (Flow % wall_dist % prec)
  call Control_Mod_Tolerance_For_Wall_Distance_Solver       &
                                                     (Flow % wall_dist % tol)
  call Control_Mod_Max_Iterations_For_Wall_Distance_Solver  &
                                                     (Flow % wall_dist % mniter)

  !--------------------------!
  !   Related to potential   !  (for flow field initialization)
  !--------------------------!
  call Control_Mod_Solver_For_Potential               (Flow % pot % solver)
  call Control_Mod_Preconditioner_For_System_Matrix   (Flow % pot % prec)
  call Control_Mod_Tolerance_For_Potential_Solver     (Flow % pot % tol)
  call Control_Mod_Max_Iterations_For_Potential_Solver(Flow % pot % mniter)

  !------------------------------!
  !   Related to heat transfer   !
  !------------------------------!
  if(Flow % heat_transfer) then
    call Control_Mod_Solver_For_Energy                (Flow % t % solver)
    call Control_Mod_Preconditioner_For_System_Matrix (Flow % t % prec)
    call Control_Mod_Tolerance_For_Energy_Solver      (Flow % t % tol)
    call Control_Mod_Max_Iterations_For_Energy_Solver (Flow % t % mniter)
  end if

  !--------------------------------!
  !   Related to multiphase Flow   !
  !--------------------------------!
  if(Flow % with_interface) then
    call Control_Mod_Solver_For_Vof                  (Vof % fun % solver)
    call Control_Mod_Preconditioner_For_System_Matrix(Vof % fun % prec)
    call Control_Mod_Tolerance_For_Vof_Solver        (Vof % fun % tol)
    call Control_Mod_Max_Iterations_For_Vof_Solver   (Vof % fun % mniter)
  end if

  !--------------------------------!
  !   Related to passive scalars   !
  !--------------------------------!
  do sc = 1, Flow % n_scalars
    phi => Flow % scalar(sc)
    call Control_Mod_Solver_For_Scalars               (phi % solver)
    call Control_Mod_Preconditioner_For_System_Matrix (phi % prec)
    call Control_Mod_Tolerance_For_Scalars_Solver     (phi % tol)
    call Control_Mod_Max_Iterations_For_Scalars_Solver(phi % mniter)
  end do

  !------------------------------!
  !   All turbuelnt quantities   !
  !------------------------------!
  do i = 1, 12
    if(i .eq.  1) tq => turb % kin
    if(i .eq.  2) tq => turb % eps
    if(i .eq.  3) tq => turb % zeta
    if(i .eq.  4) tq => turb % f22
    if(i .eq.  5) tq => turb % vis
    if(i .eq.  6) tq => turb % t2
    if(i .eq.  7) tq => turb % uu
    if(i .eq.  8) tq => turb % vv
    if(i .eq.  9) tq => turb % ww
    if(i .eq. 10) tq => turb % uv
    if(i .eq. 11) tq => turb % uw
    if(i .eq. 12) tq => turb % vw
    call Control_Mod_Solver_For_Turbulence               (tq % solver)
    call Control_Mod_Preconditioner_For_System_Matrix    (tq % prec)
    call Control_Mod_Tolerance_For_Turbulence_Solver     (tq % tol)
    call Control_Mod_Max_Iterations_For_Turbulence_Solver(tq % mniter)
  end do

  end subroutine

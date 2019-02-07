!==============================================================================!
  subroutine Read_Control_Numerical(flow)
!------------------------------------------------------------------------------!
!   Reads details about numerical models from control file.                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod,      only: Field_Type, heat_transfer
  use Var_Mod,        only: Var_Type
  use Rans_Mod,       only: kin, eps, zeta, f22, vis, t2
  use Turbulence_Mod, only: uu, vv, ww, uv, uw, vw
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
!----------------------------------[Locals]------------------------------------!
  type(Var_Type),   pointer :: pp, t, tq, ui
  character(len=80)         :: name
  integer                   :: i
!==============================================================================!

  ! Take aliases
  pp => flow % pp
  t  => flow % t

  !-------------------------!
  !   Related to momentum   !
  !-------------------------!
  do i = 1, 3
    if(i .eq. 1) ui => flow % u
    if(i .eq. 2) ui => flow % v
    if(i .eq. 3) ui => flow % w
    ui % urf   = 0.8
    ui % niter = 5
    call Control_Mod_Advection_Scheme_For_Momentum      (ui % adv_scheme)
    call Control_Mod_Blending_Coefficient_For_Momentum  (ui % blend)
    call Control_Mod_Time_Integration_Scheme            (ui % td_scheme)
    call Control_Mod_Simple_Underrelaxation_For_Momentum(ui % urf)
    call Control_Mod_Tolerance_For_Momentum_Solver      (ui % tol)
    call Control_Mod_Preconditioner_For_System_Matrix   (ui % precond)
    call Control_Mod_Max_Iterations_For_Momentum_Solver (ui % niter)
  end do

  !-------------------------!
  !   Related to pressure   !
  !-------------------------!
  pp % niter = 40
  call Control_Mod_Tolerance_For_Pressure_Solver      (pp % tol)
  call Control_Mod_Preconditioner_For_System_Matrix   (pp % precond)
  call Control_Mod_Max_Iterations_For_Pressure_Solver (pp % niter)
  call Control_Mod_Simple_Underrelaxation_For_Pressure(pp % urf)

  !------------------------------!
  !   Related to heat transfer   !
  !------------------------------!
  if(heat_transfer) then
    t % urf   = 0.7
    t % niter = 5
    call Control_Mod_Advection_Scheme_For_Energy      (t % adv_scheme)
    call Control_Mod_Blending_Coefficient_For_Energy  (t % blend)
    call Control_Mod_Time_Integration_Scheme          (t % td_scheme)
    call Control_Mod_Simple_Underrelaxation_For_Energy(t % urf)
    call Control_Mod_Tolerance_For_Energy_Solver      (t % tol)
    call Control_Mod_Preconditioner_For_System_Matrix (t % precond)
    call Control_Mod_Max_Iterations_For_Energy_Solver (t % niter)
  end if

  !------------------------------!
  !   All turbuelnt quantities   !
  !------------------------------!
  do i = 1, 12
    if(i .eq.  1) tq => kin
    if(i .eq.  2) tq => eps
    if(i .eq.  3) tq => zeta
    if(i .eq.  4) tq => f22
    if(i .eq.  5) tq => vis
    if(i .eq.  6) tq => t2
    if(i .eq.  7) tq => uu
    if(i .eq.  8) tq => vv
    if(i .eq.  9) tq => ww
    if(i .eq. 10) tq => uv
    if(i .eq. 11) tq => uw
    if(i .eq. 12) tq => vw
    tq % urf   = 1.0
    tq % niter = 6
    call Control_Mod_Advection_Scheme_For_Turbulence      (tq % adv_scheme)
    call Control_Mod_Blending_Coefficient_For_Turbulence  (tq % blend)
    call Control_Mod_Time_Integration_Scheme              (tq % td_scheme)
    call Control_Mod_Simple_Underrelaxation_For_Turbulence(tq % urf)
    call Control_Mod_Tolerance_For_Turbulence_Solver      (tq % tol)
    call Control_Mod_Preconditioner_For_System_Matrix     (tq % precond)
    call Control_Mod_Max_Iterations_For_Turbulence_Solver (tq % niter)
  end do

  end subroutine

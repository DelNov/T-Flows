!==============================================================================!
  subroutine Read_Control_Numerical(flow, turb, mult)
!------------------------------------------------------------------------------!
!   Reads details about numerical models from control file.                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod,      only: Field_Type, heat_transfer
  use Var_Mod,        only: Var_Type
  use Turb_Mod,       only: Turb_Type
  use Multiphase_Mod, only: Multiphase_Type, multiphase_model,   &
                            VOLUME_OF_FLUID
  use Control_Mod
  use Numerics_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
!----------------------------------[Locals]------------------------------------!
  type(Var_Type),   pointer :: tq, ui, phi
  character(len=80)         :: name
  integer                   :: i, sc
!==============================================================================!

  !-------------------------!
  !   Related to momentum   !
  !-------------------------!
  do i = 1, 3
    if(i .eq. 1) ui => flow % u
    if(i .eq. 2) ui => flow % v
    if(i .eq. 3) ui => flow % w
    ui % urf   = 0.8
    ui % niter = 5
    call Control_Mod_Advection_Scheme_For_Momentum            (name)
    ui % adv_scheme = Numerics_Mod_Advection_Scheme_Code      (name)
    call Control_Mod_Time_Integration_Scheme                  (name)
    ui % td_scheme = Numerics_Mod_Time_Integration_Scheme_Code(name)
    call Control_Mod_Blending_Coefficient_For_Momentum        (ui % blend)
    call Control_Mod_Simple_Underrelaxation_For_Momentum      (ui % urf)
    call Control_Mod_Tolerance_For_Momentum_Solver            (ui % tol)
    call Control_Mod_Preconditioner_For_System_Matrix         (ui % precond)
    call Control_Mod_Max_Iterations_For_Momentum_Solver       (ui % niter)
  end do

  !-------------------------!
  !   Related to pressure   !
  !-------------------------!
  flow % pp % niter = 40
  call Control_Mod_Tolerance_For_Pressure_Solver      (flow % pp % tol)
  call Control_Mod_Preconditioner_For_System_Matrix   (flow % pp % precond)
  call Control_Mod_Max_Iterations_For_Pressure_Solver (flow % pp % niter)
  call Control_Mod_Simple_Underrelaxation_For_Pressure(flow % pp % urf)

  !------------------------------!
  !   Related to heat transfer   !
  !------------------------------!
  if(heat_transfer) then
    flow % t % urf   = 0.7
    flow % t % niter = 5
    call Control_Mod_Advection_Scheme_For_Energy                    (name)
    flow % t % adv_scheme = Numerics_Mod_Advection_Scheme_Code      (name)
    call Control_Mod_Time_Integration_Scheme                        (name)
    flow % t % td_scheme = Numerics_Mod_Time_Integration_Scheme_Code(name)
    call Control_Mod_Blending_Coefficient_For_Energy  (flow % t % blend)
    call Control_Mod_Simple_Underrelaxation_For_Energy(flow % t % urf)
    call Control_Mod_Tolerance_For_Energy_Solver      (flow % t % tol)
    call Control_Mod_Preconditioner_For_System_Matrix (flow % t % precond)
    call Control_Mod_Max_Iterations_For_Energy_Solver (flow % t % niter)
  end if

  !--------------------------------!
  !   Related to multiphase flow   !
  !--------------------------------!
  if(multiphase_model .eq. VOLUME_OF_FLUID) then
    mult % vof % urf   = 0.7
    mult % vof % niter = 5
    call Control_Mod_Advection_Scheme_For_Multiphase                  (name)
    mult % vof % adv_scheme = Numerics_Mod_Advection_Scheme_Code      (name)
    call Control_Mod_Time_Integration_Scheme                          (name)
    mult % vof % td_scheme = Numerics_Mod_Time_Integration_Scheme_Code(name)
    call Control_Mod_Blending_Coefficient_For_Multiphase  (mult % vof % blend)
    call Control_Mod_Simple_Underrelaxation_For_Multiphase(mult % vof % urf)
    call Control_Mod_Tolerance_For_Multiphase_Solver      (mult % vof % tol)
    call Control_Mod_Preconditioner_For_System_Matrix     (mult % vof % precond)
    call Control_Mod_Max_Iterations_For_Multiphase_Solver (mult % vof % niter)
    ! Max Courant number and Max substep cycles
    call Control_Mod_Max_Courant_Vof(mult % courant_max_param)
    call Control_Mod_Max_Substep_Cycles_Vof(mult % n_sub_param)
    ! Max correction cycles for beta
    call Control_Mod_Max_Correction_Cycles_Beta_Vof(mult % corr_num_max)
    ! Max number convolution/smoothing steps for curvature and normal
    call Control_Mod_Max_Smoothing_Cycles_Curvature_Vof(mult % n_conv_curv)
    call Control_Mod_Max_Smoothing_Cycles_Normal_Vof(mult % n_conv_norm)
    ! Parameters for distance function
    call Control_Mod_Factor_Fictitious_Time_Vof(mult % c_tau)
    call Control_Mod_Factor_Number_Cells_Distance_Function_Vof(mult % c_eps)
    call Control_Mod_Distance_Function_Time_Integration_Scheme(name)
    mult % t_dist_scheme = Numerics_Mod_Time_Integration_Scheme_Code(name)
  end if

  !--------------------------------!
  !   Related to passive scalars   !
  !--------------------------------!
  do sc = 1, flow % n_scalars
    phi => flow % scalar(sc)
    phi % urf   = 0.7
    phi % niter = 5
    call Control_Mod_Advection_Scheme_For_Scalars              (name)
    phi % adv_scheme = Numerics_Mod_Advection_Scheme_Code      (name)
    call Control_Mod_Time_Integration_Scheme                   (name)
    phi % td_scheme = Numerics_Mod_Time_Integration_Scheme_Code(name)
    call Control_Mod_Blending_Coefficient_For_Scalars          (phi % blend)
    call Control_Mod_Simple_Underrelaxation_For_Scalars        (phi % urf)
    call Control_Mod_Tolerance_For_Scalars_Solver              (phi % tol)
    call Control_Mod_Preconditioner_For_System_Matrix          (phi % precond)
    call Control_Mod_Max_Iterations_For_Scalars_Solver         (phi % niter)
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
    tq % urf   = 1.0
    tq % niter = 6
    call Control_Mod_Advection_Scheme_For_Turbulence          (name)
    tq % adv_scheme = Numerics_Mod_Advection_Scheme_Code      (name)
    call Control_Mod_Time_Integration_Scheme                  (name)
    tq % td_scheme = Numerics_Mod_Time_Integration_Scheme_Code(name)
    call Control_Mod_Blending_Coefficient_For_Turbulence      (tq % blend)
    call Control_Mod_Simple_Underrelaxation_For_Turbulence    (tq % urf)
    call Control_Mod_Tolerance_For_Turbulence_Solver          (tq % tol)
    call Control_Mod_Preconditioner_For_System_Matrix         (tq % precond)
    call Control_Mod_Max_Iterations_For_Turbulence_Solver     (tq % niter)
  end do

  end subroutine

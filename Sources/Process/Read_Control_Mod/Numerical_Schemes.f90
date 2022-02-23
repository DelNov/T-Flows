!==============================================================================!
  subroutine Numerical_Schemes(Rc, Flow, turb, Vof, Sol)
!------------------------------------------------------------------------------!
!   Reads details about numerical schemes from control file.                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Read_Control_Type)  :: Rc
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: turb
  type(Vof_Type),    target :: Vof
  type(Solver_Type), target :: Sol
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: tq, ui, phi
  character(SL)            :: name
  integer                  :: i, sc
!==============================================================================!

  ! Take alias
  Grid => Flow % pnt_grid

  !-----------------------------------------------------!
  !   Linear solvers you want to use; native or PETSc   !
  !-----------------------------------------------------!
  call Control_Mod_Linear_Solvers(name, .true.)
  Sol % solvers = Numerics_Mod_Linear_Solvers_Code(name)

  !------------------------------------------!
  !   Pressure velocity coupling algorithm   !
  !------------------------------------------!

  ! Basic algorithm for pressure-velocity coupling (SIMPLE, PISO)
  call Control_Mod_Pressure_Momentum_Coupling(name, .true.)
  Flow % p_m_coupling = Numerics_Mod_Pressure_Momentum_Coupling_Code(name)

  if( Flow % p_m_coupling .eq. PISO) then
    call Control_Mod_Number_Of_Piso_Corrections(Flow % n_piso_corrections)
    Flow % piso_status = .false.
  end if

  ! Improvements to Rhie and Chow method (Choi, Gu)
  call Control_Mod_Choi_Correction(Flow % choi_correction, .false.)
  call Control_Mod_Gu_Correction  (Flow % gu_correction,   .false.)

  !----------------------------------!
  !   Gradient computation methods   !
  !----------------------------------!

  ! Tolerance and max iterations for computation of gradients with Gauss method
  call Control_Mod_Tolerance_For_Gauss_Gradients(Flow % gauss_tol, .false.)
  Flow % gauss_miter = 10
  call Control_Mod_Max_Gauss_Gradients_Iterations(Flow % gauss_miter,  &
                                                  .false.)
  Flow % least_miter =  4
  call Control_Mod_Max_Least_Squares_Gradients_Iterations(Flow % least_miter,  &
                                                          .false.)

  !-------------------------!
  !   Related to momentum   !
  !-------------------------!
  do i = 1, 3
    if(i .eq. 1) ui => Flow % u
    if(i .eq. 2) ui => Flow % v
    if(i .eq. 3) ui => Flow % w
    ui % urf    = 0.8
    ui % mniter = 5
    call Control_Mod_Advection_Scheme_For_Momentum            (name)
    ui % adv_scheme = Numerics_Mod_Advection_Scheme_Code      (name)
    call Control_Mod_Time_Integration_Scheme                  (name)
    ui % td_scheme = Numerics_Mod_Time_Integration_Scheme_Code(name)
    call Control_Mod_Blending_Coefficient_For_Momentum  (ui % blend)
    call Control_Mod_Simple_Underrelaxation_For_Momentum(ui % urf)
    ui % solver = 'bicg'
    call Control_Mod_Preconditioner_For_System_Matrix   (ui % prec)
    call Control_Mod_Tolerance_For_Momentum_Solver      (ui % tol)
    call Control_Mod_Max_Iterations_For_Momentum_Solver (ui % mniter)
    call Control_Mod_Gradient_Method_For_Momentum       (name)
    ui % grad_method = Numerics_Mod_Gradient_Method_Code(name)
  end do

  !-------------------------!
  !   Related to pressure   !
  !-------------------------!
  Flow % pp % mniter = 60
  Flow % pp % solver = 'cg'
  call Control_Mod_Preconditioner_For_System_Matrix   (Flow % pp % prec)
  call Control_Mod_Tolerance_For_Pressure_Solver      (Flow % pp % tol)
  call Control_Mod_Max_Iterations_For_Pressure_Solver (Flow % pp % mniter)
  call Control_Mod_Simple_Underrelaxation_For_Pressure(Flow % pp % urf)
  call Control_Mod_Gradient_Method_For_Pressure              (name)
  Flow % p  % grad_method = Numerics_Mod_Gradient_Method_Code(name)
  Flow % pp % grad_method = Numerics_Mod_Gradient_Method_Code(name)

  !------------------------------!
  !   Related to wall distance   !
  !------------------------------!
  Flow % wall_dist % mniter = 99
  Flow % wall_dist % solver = 'bicg'
  Flow % wall_dist % tol    =  1.0e-5
  call Control_Mod_Preconditioner_For_System_Matrix         &
       (Flow % wall_dist % prec)
  call Control_Mod_Tolerance_For_Wall_Distance_Solver       &
       (Flow % wall_dist % tol)
  call Control_Mod_Max_Iterations_For_Wall_Distance_Solver  &
       (Flow % wall_dist % mniter)
  call Control_Mod_Gradient_Method_For_Wall_Distance(name)

  !--------------------------!
  !   Related to potential   !  (for flow field initialization)
  !--------------------------!
  Flow % pot % mniter = 99
  Flow % pot % solver = 'bicg'
  Flow % pot % tol    =  1.0e-5
  call Control_Mod_Preconditioner_For_System_Matrix   (Flow % pot % prec)
  call Control_Mod_Tolerance_For_Potential_Solver     (Flow % pot % tol)
  call Control_Mod_Max_Iterations_For_Potential_Solver(Flow % pot % mniter)

  !------------------------------!
  !   Related to heat transfer   !
  !------------------------------!
  if(Flow % heat_transfer) then
    Flow % t % urf    = 0.7
    Flow % t % mniter = 5
    call Control_Mod_Advection_Scheme_For_Energy                    (name)
    Flow % t % adv_scheme = Numerics_Mod_Advection_Scheme_Code      (name)
    call Control_Mod_Time_Integration_Scheme                        (name)
    Flow % t % td_scheme = Numerics_Mod_Time_Integration_Scheme_Code(name)
    call Control_Mod_Blending_Coefficient_For_Energy  (Flow % t % blend)
    call Control_Mod_Simple_Underrelaxation_For_Energy(Flow % t % urf)
    Flow % t % solver = 'bicg'
    call Control_Mod_Preconditioner_For_System_Matrix (Flow % t % prec)
    call Control_Mod_Tolerance_For_Energy_Solver      (Flow % t % tol)
    call Control_Mod_Max_Iterations_For_Energy_Solver (Flow % t % mniter)
    call Control_Mod_Gradient_Method_For_Energy               (name)
    Flow % t % grad_method = Numerics_Mod_Gradient_Method_Code(name)
  end if

  !--------------------------------!
  !   Related to multiphase Flow   !
  !--------------------------------!
  if(Flow % with_interface) then
    Vof % fun % urf    = 0.7
    Vof % fun % mniter = 5
    call Control_Mod_Advection_Scheme_For_Vof                        (name)
    Vof % fun % adv_scheme = Numerics_Mod_Advection_Scheme_Code      (name)
    call Control_Mod_Time_Integration_Scheme                         (name)
    Vof % fun % td_scheme = Numerics_Mod_Time_Integration_Scheme_Code(name)
    call Control_Mod_Blending_Coefficient_For_Vof    (Vof % fun % blend)
    call Control_Mod_Simple_Underrelaxation_For_Vof  (Vof % fun % urf)
    Vof % fun % solver = 'bicg'
    call Control_Mod_Preconditioner_For_System_Matrix(Vof % fun % prec)
    call Control_Mod_Tolerance_For_Vof_Solver        (Vof % fun % tol)
    call Control_Mod_Max_Iterations_For_Vof_Solver   (Vof % fun % mniter)
    ! Max Courant number and Max substep cycles
    call Control_Mod_Max_Courant_Vof       (Vof % courant_max_param)
    call Control_Mod_Max_Substep_Cycles_Vof(Vof % n_sub_param)
    call Control_Mod_Max_Correction_Cycles_Beta_Vof    (Vof % corr_num_max)
    call Control_Mod_Max_Smoothing_Cycles_Curvature_Vof(Vof % n_conv_curv)
    call Control_Mod_Max_Smoothing_Cycles_Normal_Vof   (Vof % n_conv_norm)
    call Control_Mod_Skewness_Correction_Vof           (Vof % skew_corr)
    call Control_Mod_Gradient_Method_For_Vof                   (name)
    Vof % fun % grad_method = Numerics_Mod_Gradient_Method_Code(name)
  end if

  !--------------------------------!
  !   Related to passive scalars   !
  !--------------------------------!
  do sc = 1, Flow % n_scalars
    phi => Flow % scalar(sc)
    phi % urf    = 0.7
    phi % mniter = 5
    call Control_Mod_Advection_Scheme_For_Scalars              (name)
    phi % adv_scheme = Numerics_Mod_Advection_Scheme_Code      (name)
    call Control_Mod_Time_Integration_Scheme                   (name)
    phi % td_scheme = Numerics_Mod_Time_Integration_Scheme_Code(name)
    call Control_Mod_Blending_Coefficient_For_Scalars  (phi % blend)
    call Control_Mod_Simple_Underrelaxation_For_Scalars(phi % urf)
    phi % solver = 'bicg'
    call Control_Mod_Preconditioner_For_System_Matrix  (phi % prec)
    call Control_Mod_Tolerance_For_Scalars_Solver      (phi % tol)
    call Control_Mod_Max_Iterations_For_Scalars_Solver (phi % mniter)
    call Control_Mod_Gradient_Method_For_Scalars         (name)
    phi % grad_method = Numerics_Mod_Gradient_Method_Code(name)
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
    tq % urf    = 1.0
    tq % mniter = 6
    call Control_Mod_Advection_Scheme_For_Turbulence          (name)
    tq % adv_scheme = Numerics_Mod_Advection_Scheme_Code      (name)
    call Control_Mod_Time_Integration_Scheme                  (name)
    tq % td_scheme = Numerics_Mod_Time_Integration_Scheme_Code(name)
    call Control_Mod_Blending_Coefficient_For_Turbulence  (tq % blend)
    call Control_Mod_Simple_Underrelaxation_For_Turbulence(tq % urf)
    tq % solver = 'bicg'
    call Control_Mod_Preconditioner_For_System_Matrix    (tq % prec)
    call Control_Mod_Tolerance_For_Turbulence_Solver     (tq % tol)
    call Control_Mod_Max_Iterations_For_Turbulence_Solver(tq % mniter)
    call Control_Mod_Gradient_Method_For_Turbulence     (name)
    tq % grad_method = Numerics_Mod_Gradient_Method_Code(name)
  end do

  end subroutine

!==============================================================================!
  subroutine Numerical_Schemes(Rc, Flow, Turb, Vof)
!------------------------------------------------------------------------------!
!   Reads details about numerical schemes from control file.                   !
!                                                                              !
!   Good practice: default values, outlined in Documents/all_control_keywords, !
!   should be defined only in Control_Mod, not be scattered around the code.   !
!   In other words, Control_Mod changes less frequently than other parts of    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Read_Controls_Type), intent(in) :: Rc
  type(Field_Type), target              :: Flow
  type(Turb_Type),  target              :: Turb
  type(Vof_Type),   target              :: Vof
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: tq, ui, phi
  character(SL)            :: name
  integer                  :: i, sc
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Rc)
!==============================================================================!

  ! Take alias
  Grid => Flow % pnt_grid

  !------------------------------------------!
  !   Pressure velocity coupling algorithm   !
  !------------------------------------------!

  ! Basic algorithm for pressure-velocity coupling (SIMPLE, PISO)
  call Control % Pressure_Momentum_Coupling(name, .true.)
  Flow % p_m_coupling = Numerics_Mod_Pressure_Momentum_Coupling_Code(name)

  if( Flow % p_m_coupling .eq. PISO) then
    call Control % Number_Of_Piso_Corrections(Flow % n_piso_corrections)
    Flow % inside_piso_loop = .false.
  end if

  ! Improvements to Rhie and Chow method (Choi, Gu)
  call Control % Choi_Correction(Flow % choi_correction, .false.)
  call Control % Gu_Correction  (Flow % gu_correction,   .false.)

  ! Report volume balance (in a separate file)
  call Control % Report_Volume_Balance(Flow % rep_vol_balance, .false.)

  !----------------------------------!
  !   Gradient computation methods   !
  !----------------------------------!

  ! Tolerance and max iterations for computation of gradients with Gauss method
  call Control % Tolerance_For_Gauss_Gradients (Flow % gauss_tol,   .false.)
  call Control % Max_Gauss_Gradients_Iterations(Flow % gauss_miter, .false.)
  call Control % Max_Least_Squares_Gradients_Iterations(Flow % least_miter,  &
                                                          .false.)
  !-------------------------!
  !   Related to momentum   !
  !-------------------------!
  do i = 1, 3
    if(i .eq. 1) ui => Flow % u
    if(i .eq. 2) ui => Flow % v
    if(i .eq. 3) ui => Flow % w
    call Control % Advection_Scheme_For_Momentum              (name)
    ui % adv_scheme = Numerics_Mod_Advection_Scheme_Code      (name)
    call Control % Time_Integration_Scheme                    (name)
    ui % td_scheme = Numerics_Mod_Time_Integration_Scheme_Code(name)
    call Control % Blending_Coefficient_For_Momentum    (ui % blend)
    call Control % Simple_Underrelaxation_For_Momentum  (ui % urf)
    call Control % Gradient_Method_For_Momentum         (name)
    ui % grad_method = Numerics_Mod_Gradient_Method_Code(name)
    call Control % Blend_System_Matrices(ui % blend_matrix, .false.)
  end do

  !-------------------------!
  !   Related to pressure   !
  !-------------------------!
  call Control % Simple_Underrelaxation_For_Pressure(Flow % pp % urf)
  call Control % Gradient_Method_For_Pressure                (name)
  Flow % p  % grad_method = Numerics_Mod_Gradient_Method_Code(name)
  Flow % pp % grad_method = Numerics_Mod_Gradient_Method_Code(name)
  call Control % Blend_System_Matrices(Flow % pp % blend_matrix, .false.)

  !------------------------------!
  !   Related to wall distance   !
  !------------------------------!
  call Control % Gradient_Method_For_Wall_Distance                  (name)
  Flow % wall_dist % grad_method = Numerics_Mod_Gradient_Method_Code(name)
  call Control % Blend_System_Matrices(Flow % wall_dist % blend_matrix, .false.)

  !--------------------------!
  !   Related to potential   !  (for flow field initialization, nothing here)
  !--------------------------!


  !------------------------------!
  !   Related to heat transfer   !
  !------------------------------!
  if(Flow % heat_transfer) then
    call Control % Advection_Scheme_For_Energy                      (name)
    Flow % t % adv_scheme = Numerics_Mod_Advection_Scheme_Code      (name)
    call Control % Time_Integration_Scheme                          (name)
    Flow % t % td_scheme = Numerics_Mod_Time_Integration_Scheme_Code(name)
    call Control % Blending_Coefficient_For_Energy  (Flow % t % blend)
    call Control % Simple_Underrelaxation_For_Energy(Flow % t % urf)
    call Control % Gradient_Method_For_Energy                 (name)
    Flow % t % grad_method = Numerics_Mod_Gradient_Method_Code(name)
    call Control % Blend_System_Matrices(Flow % t % blend_matrix, .false.)
  end if

  !--------------------------------!
  !   Related to multiphase Flow   !
  !--------------------------------!
  if(Flow % with_interface) then
    call Control % Advection_Scheme_For_Vof                          (name)
    Vof % fun % adv_scheme = Numerics_Mod_Advection_Scheme_Code      (name)
    call Control % Time_Integration_Scheme                           (name)
    Vof % fun % td_scheme = Numerics_Mod_Time_Integration_Scheme_Code(name)
    call Control % Blending_Coefficient_For_Vof  (Vof % fun % blend)
    call Control % Simple_Underrelaxation_For_Vof(Vof % fun % urf)
    ! Max Courant number and Max substep cycles
    call Control % Max_Courant_Vof                   (Vof % courant_max_param)
    call Control % Max_Substep_Cycles_Vof            (Vof % n_sub_param)
    call Control % Max_Correction_Cycles_Beta_Vof    (Vof % corr_num_max)
    call Control % Max_Smoothing_Cycles_Curvature_Vof(Vof % n_conv_curv)
    call Control % Max_Smoothing_Cycles_Normal_Vof   (Vof % n_conv_norm)
    call Control % Skewness_Correction_Vof           (Vof % skew_corr)
    call Control % Gradient_Method_For_Vof                     (name)
    Vof % fun % grad_method = Numerics_Mod_Gradient_Method_Code(name)
    call Control % Blend_System_Matrices(Vof % fun % blend_matrix, .false.)
  end if

  !--------------------------------!
  !   Related to passive scalars   !
  !--------------------------------!
  do sc = 1, Flow % n_scalars
    phi => Flow % scalar(sc)
    call Control % Advection_Scheme_For_Scalars                (name)
    phi % adv_scheme = Numerics_Mod_Advection_Scheme_Code      (name)
    call Control % Time_Integration_Scheme                     (name)
    phi % td_scheme = Numerics_Mod_Time_Integration_Scheme_Code(name)
    call Control % Blending_Coefficient_For_Scalars            (phi % blend)
    call Control % Simple_Underrelaxation_For_Scalars          (phi % urf)
    call Control % Gradient_Method_For_Scalars                 (name)
    phi % grad_method = Numerics_Mod_Gradient_Method_Code      (name)
    call Control % Blend_System_Matrices(phi % blend_matrix, .false.)
  end do

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
    call Control % Advection_Scheme_For_Turbulence            (name)
    tq % adv_scheme = Numerics_Mod_Advection_Scheme_Code      (name)
    call Control % Time_Integration_Scheme                    (name)
    tq % td_scheme = Numerics_Mod_Time_Integration_Scheme_Code(name)
    call Control % Blending_Coefficient_For_Turbulence        (tq % blend)
    call Control % Simple_Underrelaxation_For_Turbulence      (tq % urf)
    call Control % Gradient_Method_For_Turbulence             (name)
    tq % grad_method = Numerics_Mod_Gradient_Method_Code      (name)
    call Control % Blend_System_Matrices(tq % blend_matrix, .false.)
  end do

  end subroutine

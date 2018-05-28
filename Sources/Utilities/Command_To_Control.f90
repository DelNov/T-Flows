include '../Shared/To_Lower_Case.f90'
include '../Shared/To_Upper_Case.f90'
include 'Tokenizer.f90'

!==============================================================================!
  program Command_To_Control
!------------------------------------------------------------------------------!
!   Reads command file and turns it into control file.                         !
!------------------------------------------------------------------------------!
  use Tokenizer_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: answer_name
  character(len=80) :: name_in
  character(len=80) :: answer, model
  integer           :: it, i, n_mon
  logical           :: HOT, SIMPLE
!==============================================================================!

  HOT    = .false.
  SIMPLE = .false.

  !---------------------------!
  !   Open the command file   !
  !---------------------------!
  open(CMN_FILE, file='T-FlowS.cmn')
  cmn_line_count = 0

  !-----------------------!
  !   Read answer name   !
  !-----------------------!
  print *, '# Input answer name:'
  call Tokenizer_Mod_Read_Line(CMN_FILE)

  ! Print entry for the control file
  print *, 'PROBLEM_NAME    ', trim(line % tokens(1))

  !-------------------------------!
  !   Reads type of the answer   !
  !-------------------------------!
  ! call Read_Problem           ! bad practice, should be avoided
  print *, '# Type of answer: '
  print *, '# CHANNEL          -> Channel flow'
  print *, '# PIPE             -> Pipe flow'
  print *, '# JET              -> Impinging jet flow'
  print *, '# TEST             -> Test Laplacian equation'
  print *, '# OTHER            -> All the other answers'
  print *, '# ROUGH            -> Problems with roughness'
  print *, '# HOT              -> Problems with temperature'
  print *, '# XHOM, YHOM, ZHOM -> Homogeneous directions'
  print *, '# TGV -> Taylor-Green Vortex test case'
  print *, '# BUOY -> Buoyancy flows (Automatically turns HOT on)'
  print *, '# RB_CONV          -> Rayleigh-Barnard convection'
  print *, '# URANS            -> Unsteady RANS'
  print *, '# BACKSTEP         -> Backstep flow'
  call Tokenizer_Mod_Read_Line(CMN_FILE)

  do it = 1, line % n_tokens
    read(line % tokens(it),'(A8)')  answer
    call To_Lower_Case(answer)

    ! Print entry for the control file
    if(it .eq. 1) print *, 'PROBLEM_TYPE    ', trim(answer)

    if(answer .eq. 'hot') then

      ! Print entry for the control file
      print *, 'HEAT_TRANSFER    yes'
      HOT = .true.

    else if(answer .eq. 'urans') then

      ! Print entry for the control file
      print *, 'TURBULENCE_MODEL_VARIANT    urans'

    else if(answer .eq. 'rot') then
      print *, '# Angular velocity vector: '
      call Tokenizer_Mod_Read_Line(CMN_FILE)

      ! Print entry for the control file
      print *, 'ANGULAR_VELOCITY_VECTOR    ',  trim(line % tokens(1)), '  ',  &
                                               trim(line % tokens(2)), '  ',  &
                                               trim(line % tokens(3))

    else if(answer .eq. 'rb_conv') then
      print *, '# Gravitational constant in x, y and z directions: '
      call Tokenizer_Mod_Read_Line(CMN_FILE)

      print *, 'BUOYANCY    yes'

      ! Print entry for the control file
      print *, 'GRAVITATIONAL_VECTOR    ',  trim(line % tokens(1)), '  ', &
                                            trim(line % tokens(2)), '  ', &
                                            trim(line % tokens(3))
      print *, 'REFERENCE_TEMPERATURE    ', trim(line % tokens(4))

    else if(answer .eq. 'buoy') then
      print *, '# Gravitational constant in x, y and z directions: '
      call Tokenizer_Mod_Read_Line(CMN_FILE)

      print *, 'BUOYANCY    yes'

      ! Print entry for the control file
      print *, 'GRAVITATIONAL_VECTOR    ',  trim(line % tokens(1)), '  ', &
                                            trim(line % tokens(2)), '  ', &
                                            trim(line % tokens(3))
      print *, 'REFERENCE_TEMPERATURE    ', trim(line % tokens(4))

    else if(answer .eq. 'rough') then
      print *, '# Reading roughness coefficient Zo'
      call Tokenizer_Mod_Read_Line(CMN_FILE)

      ! Print entry for the control file
      print *, 'ROUGHNESS_COEFFICIENT    ', trim(line % tokens(1))

!   else
!     print *, '# answer = ', answer
!     print *, '# Error in input ! Exiting'
!     stop
    endif
  end do

  !------------------------!
  !   Reads restart file   !
  !------------------------!
  ! call Load_Restart(grid, restar)
  print *, '# Input restart file name [skip cancels]:'
  call Tokenizer_Mod_Read_Line(CMN_FILE)
  read(line % tokens(1), '(A80)') answer
  call To_Upper_Case(answer)

  ! Print entry for the control file
  if(answer .ne. 'SKIP') then
    print *, 'LOAD_RESTART_NAME    ', trim(line % tokens(1))
  end if

  !-----------------------!
  !   Reads T-Flows.cmn   !
  !-----------------------!
  ! call ReaCom(grid, restar)
  ! The number of time steps
  print *, '#==============================================='
  print *, '# Enter the number of time steps: '
  print *, '# (type 0 if you just want to analyse results)'
  print *, '#-----------------------------------------------'
  call Tokenizer_Mod_Read_Line(CMN_FILE)

  ! Print entry for the control file
  print *, 'NUMBER_OF_TIME_STEPS    ', trim(line % tokens(1))

  ! Starting time step for statistics
  print *, '# Starting time step for statistics '
  call Tokenizer_Mod_Read_Line(CMN_FILE)

  ! Print entry for the control file
  print *, 'STARTING_TIME_STEP_FOR_STATISTICS    ', trim(line % tokens(1))

  ! We agreed that budgets won't be supported
  ! if(BUDG .eq. YES) then
  !   read(line % tokens(2),*) Nbudg
  ! end if

  print *, '# Number of monitoring points:'
  call Tokenizer_Mod_Read_Line(CMN_FILE)
  read(line % tokens(1), *) n_mon

  ! Print entry for the control file
  print *, 'NUMBER_OF_MONITORING_POINTS    ', n_mon

  print *, '# Enter the coordinates of monitoring point(s)'
  do i = 1, n_mon
    call Tokenizer_Mod_Read_Line(CMN_FILE)

    ! Print entry for the control file
    print '(a,i3.3,7a)', ' MONITORING_POINT_',          &
                         i,                      '  ',  &
                         trim(line % tokens(1)), '  ',  &
                         trim(line % tokens(2)), '  ',  &
                         trim(line % tokens(3))

  end do

  ! Plane for calcution of overall mass fluxes
  ! do m = 1, grid % n_materials
  print *, '# Enter the coordinates of monitoring plane: '
  call Tokenizer_Mod_Read_Line(CMN_FILE)

  ! Print entry for the control file
  print *, 'POINT_FOR_MONITORING_PLANES    ', trim(line % tokens(1)), '  ',  &
                                              trim(line % tokens(2)), '  ',  &
                                              trim(line % tokens(3))
  ! end do

  ! Turbulence model
  print *, '# Type of simulation: '
  print *, '# DNS      -> Direct Numerical Simulation'
  print *, '# LES      -> Large Eddy Simulation'
  print *, '# K_EPS    -> High Reynolds k-eps model.'
  print *, '# K_EPS_VV -> Durbin`s model.'
  print *, '# SPA_ALL  -> Spalart-Allmaras model.'
  print *, '# ZETA  -> k-eps-zeta-f model.'
  call Tokenizer_Mod_Read_Line(CMN_FILE)
  read(line % tokens(1),'(A)')  model

  ! Print entry for the control file
  call To_Lower_Case(model)

  if(model .eq. 'ZETA') then
    print *, 'TURBULENCE_MODEL    k_eps_zeta_f'

  else if(model .eq. 'HYB_ZETA') then
    print *, 'TURBULENCE_MODEL    hybrid_k_eps_zeta_f'

  else if(model .eq. 'K_EPS_VV') then
    print *, 'TURBULENCE_MODEL    k_eps_v2'

  else if(model .eq. 'HYB_PITM') then
    print *, 'TURBULENCE_MODEL    hybrid_pitm'

  else if(model .eq. 'DES_SAP') then
    print *, 'TURBULENCE_MODEL    des_spalart'

  else if(model .eq. 'SPA_ALL') then
    print *, 'TURBULENCE_MODEL    spalart_allmaras'

  else if(model .eq. 'EBM') then
    print *, 'TURBULENCE_MODEL    reynolds_stress_model'

  else if(model .eq. 'HJ') then
    print *, 'TURBULENCE_MODEL    hanjalic_jakirlic'

  ! K_EPS, LES, DES
  else
    print *, 'TURBULENCE_MODEL    ', trim(model)
  end if

  call To_Upper_Case(model)

  ! Prepare for checking variants of the model
  call To_Upper_Case(line % tokens(2))

  ! Variants of RSM
  if(model .eq. 'EBM' .or. model .eq. 'HJ') then
    if(line % tokens(2) .eq. 'HYB') then
      ! Print entry for the control file
      print *, 'TURBULENCE_MODEL_VARIANT    hybrid'
    else
      ! Print entry for the control file
      print *, 'TURBULENCE_MODEL_VARIANT    pure'
    end if
  end if

  ! Variants of K_EPS
  if(model .eq. 'K_EPS') then
    if(line % tokens(2) .eq. 'LRE') then
      ! Print entry for the control file
      print *, 'TURBULENCE_MODEL_VARIANT    low_re'
    else if(line % tokens(2) .eq. 'HRE') then
      ! Print entry for the control file
      print *, 'TURBULENCE_MODEL_VARIANT    high_re'
    end if
  end if

  if(model .eq. 'LES' .or. model .eq. 'HYB_ZETA') then
    if(line % tokens(2) .eq. 'SMAG') then
      print *, 'TURBULENCE_MODEL_VARIANT    smagorinsky'
      print *, 'TURBULENCE_MODEL_CONSTANT    ', trim(line % tokens(3))
    else if(line % tokens(2) .eq. 'DYN') then
      print *, 'TURBULENCE_MODEL_VARIANT    dynamic'
    else if(line % tokens(2) .eq. 'WALE') then
      print *, 'TURBULENCE_MODEL_VARIANT    wale'
    end if
  end if

!@do m = 1, grid % n_materials
  if(model .eq. 'LES'       .or.  &
     model .eq. 'DNS'       .or.  &
     model .eq. 'DES_SPA'   .or.  &
     model .eq. 'HYB_PITM'  .or.  &
     model .eq. 'HYB_ZETA') then
     print *, '# Do you want to shake the velocity field ?'
     print *, '# YES -> shake'
     print *, '# NO  -> do not shake'
     call Tokenizer_Mod_Read_Line(CMN_FILE)
     call To_Upper_Case(line % tokens(1))

     if(line % tokens(1) .eq. 'YES') then
       print *, 'PERTURB_MOMENTUM    yes'
       print *, '# For how many time steps you want to shake ?'
       call Tokenizer_Mod_Read_Line(CMN_FILE)
       print *, 'MOMENTUM_PERTURBATION_UNTIL    ', trim(line % tokens(1))
       print *, '# Interval for shaking:'
       call Tokenizer_Mod_Read_Line(CMN_FILE)
       print *, 'MOMENTUM_PERTURBATION_INTERVAL    ', trim(line % tokens(1))
     else
       print *, 'PERTURB_MOMENTUM    no'
     end if
  end if
!@end do

  ! Time stepping scheme
  print *, '# Algorythm for time-integration: '
  print *, '# SIMPLE [Nini] -> S. I. M. P. L. E.'
  print *, '# FRACTION      -> Fractional step method'
  call Tokenizer_Mod_Read_Line(CMN_FILE)
  call To_Upper_Case(line % tokens(1))
  if(line % tokens(1) .eq. 'FRACTION') then
    print *, 'PRESSURE_MOMENTUM_COUPLING    projection'
  else if(line % tokens(1) .eq. 'SIMPLE') then

    SIMPLE = .true.

    print *, 'PRESSURE_MOMENTUM_COUPLING    simple'
    if(line % n_tokens .eq. 2) then
      print *, 'MAX_NUMBER_OF_SIMPLE_ITERATIONS    ', trim(line % tokens(2))
    end if

    print *, '# Under Relaxation Factor for velocity'
    call Tokenizer_Mod_Read_Line(CMN_FILE)
    print *, 'SIMPLE_UNDERRELAXATION_FOR_MOMENTUM    ', trim(line % tokens(1))
    print *, '# Under Relaxation Factor for pressure'
    call Tokenizer_Mod_Read_Line(CMN_FILE)
    print *, 'SIMPLE_UNDERRELAXATION_FOR_PRESSURE    ', trim(line % tokens(1))
    if(HOT) then
      print *, '# Under Relaxation Factor for temperature'
      call Tokenizer_Mod_Read_Line(CMN_FILE)
      print *, 'SIMPLE_UNDERRELAXATION_FOR_ENERGY    ', trim(line % tokens(1))
    end if
    if(model .ne. 'LES' .and. model .ne. 'DNS') then
      print *, '# Under Relaxation Factor for turbulent variables'
      call Tokenizer_Mod_Read_Line(CMN_FILE)
      print *, 'SIMPLE_UNDERRELAXATION_FOR_TURBULENCE    ', trim(line % tokens(1))
    end if

  end if

  print *, '# Integration of inertial terms: '
  print *, '# LIN -> Linear'
  print *, '# PAR -> Parabolic'
  call Tokenizer_Mod_Read_Line(CMN_FILE)
  read(line % tokens(1),'(A)')  answer
  call To_Upper_Case(answer)
  if(answer .eq. 'LIN') then
    print *, 'TIME_INTEGRATION_FOR_INERTIA    linear'
  else if(answer .eq. 'PAR') then
    print *, 'TIME_INTEGRATION_FOR_INERTIA    parabolic'
  endif

  print *, '# Integration of advection terms: '
  print *, '# AB -> Adams-Bashforth'
  print *, '# CN -> Crank-Nicholson'
  print *, '# FI -> Fully Implicit'
  call Tokenizer_Mod_Read_Line(CMN_FILE)
  read(line % tokens(1),'(A)')  answer
  call To_Upper_Case(answer)
  if(answer .eq. 'AB') then
    print *, 'TIME_INTEGRATION_FOR_ADVECTION    adams_bashforth'
  else if(answer .eq. 'CN') then
    print *, 'TIME_INTEGRATION_FOR_ADVECTION    crank_nicolson'
  else if(answer .eq. 'FI') then
    print *, 'TIME_INTEGRATION_FOR_ADVECTION    fully_implicit'
  endif

  print *, '# Integration of diffusion terms: '
  print *, '# AB -> Adams-Bashforth'
  print *, '# CN -> Crank-Nicholson'
  print *, '# FI -> Fully Implicit'
  call Tokenizer_Mod_Read_Line(CMN_FILE)
  read(line % tokens(1),'(A)')  answer
  call To_Upper_Case(answer)
  if(answer .eq. 'AB') then
    print *, 'TIME_INTEGRATION_FOR_DIFFUSION    adams_bashforth'
  else if(answer .eq. 'CN') then
    print *, 'TIME_INTEGRATION_FOR_DIFFUSION    crank_nicolson'
  else if(answer .eq. 'FI') then
    print *, 'TIME_INTEGRATION_FOR_DIFFUSION    fully_implicit'
  endif

  print *, '# Integration of cross-diffusion terms: '
  print *, '# AB -> Adams-Bashforth'
  print *, '# CN -> Crank-Nicholson'
  print *, '# FI -> Fully Implicit'
  call Tokenizer_Mod_Read_Line(CMN_FILE)
  read(line % tokens(1),'(A)')  answer
  call To_Upper_Case(answer)
  if(answer .eq. 'AB') then
    print *, 'TIME_INTEGRATION_FOR_CROSS_DIFFUSION    adams_bashforth'
  else if(answer .eq. 'CN') then
    print *, 'TIME_INTEGRATION_FOR_CROSS_DIFFUSION    crank_nicolson'
  else if(answer .eq. 'FI') then
    print *, 'TIME_INTEGRATION_FOR_CROSS_DIFFUSION    fully_implicit'
  endif

  ! Upwind blending for momentum
  print *, '# Convetive schemes for momentum equation:'
  print *, '# Do you want to use upwind blending: '
  print *, '# YES       -> use blening'
  print *, '# NO        -> do not use blending'
  print *, '# CDS       -> central differencing'
  print *, '# LUDS      -> linear upwind'
  print *, '# QUICK     -> self descriptive'
  print *, '# MINMOD    -> self descriptive'
  print *, '# SMART     -> self descriptive'
  print *, '# AVL_SMART -> self descriptive'
  call Tokenizer_Mod_Read_Line(CMN_FILE)
  read(line % tokens(1),'(A)')  answer
  call To_Upper_Case(answer)
  if(answer .eq. 'YES') then
    if(line % n_tokens .eq. 1) then
      print *, 'ADVECTION_SCHEME_FOR_MOMENTUM    blended'
      print *, 'BLENDING_COEFFICIENT_MOMENTUM    0.5'
    else
      print *, 'ADVECTION_SCHEME_FOR_MOMENTUM    blended'
      print *, 'BLENDING_COEFFICIENT_MOMENTUM    ', trim(line % tokens(2))
    end if
  else if(answer .eq. 'NO') then
    print *, 'ADVECTION_SCHEME_FOR_MOMENTUM    central'
  else if(answer .eq. 'UDS') then
    print *, 'ADVECTION_SCHEME_FOR_MOMENTUM    upwind'
  else if(answer .eq. 'CDS') then
    print *, 'ADVECTION_SCHEME_FOR_MOMENTUM    central'
    else if(answer .eq. 'LUDS') then
    print *, 'ADVECTION_SCHEME_FOR_MOMENTUM    luds'
  else if(answer .eq. 'QUICK') then
    print *, 'ADVECTION_SCHEME_FOR_MOMENTUM    quick'
  else if(answer .eq. 'MINMOD') then
    print *, 'ADVECTION_SCHEME_FOR_MOMENTUM    minmod'
  else if(answer .eq. 'SMART') then
    print *, 'ADVECTION_SCHEME_FOR_MOMENTUM    smart'
  else if(answer .eq. 'AVL_SMART') then
    print *, 'ADVECTION_SCHEME_FOR_MOMENTUM    avl_smart'
  else if(answer .eq. 'SUPERBEE') then
    print *, 'ADVECTION_SCHEME_FOR_MOMENTUM    superbee'
  else if(answer .eq. 'GAMMA') then
    print *, 'ADVECTION_SCHEME_FOR_MOMENTUM    gamma'
  end if

  if(HOT) then
    print *, '# Convetive schemes for energy equation:'
    print *, '# Do you want to use upwind blending: '
    print *, '# YES       -> use blening'
    print *, '# NO        -> do not use blending'
    print *, '# CDS       -> central differencing'
    print *, '# LUDS      -> linear upwind'
    print *, '# QUICK     -> self descriptive'
    print *, '# MINMOD    -> self descriptive'
    print *, '# SMART     -> self descriptive'
    print *, '# AVL_SMART -> self descriptive'
    call Tokenizer_Mod_Read_Line(CMN_FILE)
    read(line % tokens(1),'(A)')  answer
    call To_Upper_Case(answer)
    if(answer .eq. 'YES') then
      if(line % n_tokens .eq. 1) then
        print *, 'ADVECTION_SCHEME_FOR_ENERGY    blended'
        print *, 'BLENDING_COEFFICIENT_ENERGY    0.5'
      else
        print *, 'ADVECTION_SCHEME_FOR_ENERGY    blended'
        print *, 'BLENDING_COEFFICIENT_ENERGY    ', trim(line % tokens(2))
      end if
    else if(answer .eq. 'NO') then
      print *, 'ADVECTION_SCHEME_FOR_ENERGY    no'  ! Whatever that means :-(
    else if(answer .eq. 'UDS') then
      print *, 'ADVECTION_SCHEME_FOR_ENERGY    upwind'
    else if(answer .eq. 'CDS') then
      print *, 'ADVECTION_SCHEME_FOR_ENERGY    central'
    else if(answer .eq. 'LUDS') then
      print *, 'ADVECTION_SCHEME_FOR_ENERGY    luds'
    else if(answer .eq. 'QUICK') then
      print *, 'ADVECTION_SCHEME_FOR_ENERGY    quick'
    else if(answer .eq. 'MINMOD') then
      print *, 'ADVECTION_SCHEME_FOR_ENERGY    minmod'
    else if(answer .eq. 'SMART') then
      print *, 'ADVECTION_SCHEME_FOR_ENERGY    smart'
    else if(answer .eq. 'AVL_SMART') then
      print *, 'ADVECTION_SCHEME_FOR_ENERGY    avl_smart'
    else if(answer .eq. 'SUPERBEE') then
      print *, 'ADVECTION_SCHEME_FOR_ENERGY    superbee'
    else if(answer .eq. 'GAMMA') then
      print *, 'ADVECTION_SCHEME_FOR_ENERGY    gamma'
    end if
  end if

  if(model .ne. 'LES' .and. model .ne. 'DNS') then
    print *, '# Convetive schemes for turbulence equation:'
    print *, '# Do you want to use upwind blending: '
    print *, '# YES       -> use blening'
    print *, '# NO        -> do not use blending'
    print *, '# CDS       -> central differencing'
    print *, '# LUDS      -> linear upwind'
    print *, '# QUICK     -> self descriptive'
    print *, '# MINMOD    -> self descriptive'
    print *, '# SMART     -> self descriptive'
    print *, '# AVL_SMART -> self descriptive'
    call Tokenizer_Mod_Read_Line(CMN_FILE)
    read(line % tokens(1),'(A)')  answer
    call To_Upper_Case(answer)
    if(answer .eq. 'YES') then
      if(line % n_tokens. eq. 1) then
        print *, 'ADVECTION_SCHEME_FOR_TURBULENCE    blended'
        print *, 'BLENDING_COEFFICIENT_TURBULENCE    0.5'
      else
        print *, 'ADVECTION_SCHEME_FOR_TURBULENCE    blended'
        print *, 'BLENDING_COEFFICIENT_TURBULENCE    ', trim(line % tokens(2))
      end if
    else if(answer .eq. 'NO') then
      print *, 'ADVECTION_SCHEME_FOR_TURBULENCE    no'  ! Whatever that means :-(
    else if(answer .eq. 'UDS') then
      print *, 'ADVECTION_SCHEME_FOR_TURBULENCE    upwind'
    else if(answer .eq. 'CDS') then
      print *, 'ADVECTION_SCHEME_FOR_TURBULENCE    central'
    else if(answer .eq. 'LUDS') then
      print *, 'ADVECTION_SCHEME_FOR_TURBULENCE    luds'
    else if(answer .eq. 'QUICK') then
      print *, 'ADVECTION_SCHEME_FOR_TURBULENCE    quick'
    else if(answer .eq. 'MINMOD') then
      print *, 'ADVECTION_SCHEME_FOR_TURBULENCE    minmod'
    else if(answer .eq. 'SMART') then
      print *, 'ADVECTION_SCHEME_FOR_TURBULENCE    smart'
    else if(answer .eq. 'AVL_SMART') then
      print *, 'ADVECTION_SCHEME_FOR_TURBULENCE    avl_smart'
    else if(answer .eq. 'SUPERBEE') then
      print *, 'ADVECTION_SCHEME_FOR_TURBULENCE    superbee'
    else if(answer .eq. 'GAMMA') then
      print *, 'ADVECTION_SCHEME_FOR_TURBULENCE    gamma'
    end if
  end if

  ! Solver parameters
  print *, '# Preconditioning of the system matrix: '
  print *, '# NO -> No preconditioning'
  print *, '# DI -> Diagonal preconditioning'
  print *, '# IC -> Incomplete Cholesky'
  call Tokenizer_Mod_Read_Line(CMN_FILE)
  read(line % tokens(1),'(A)')  answer
  call To_Upper_Case(answer)
  if(answer .eq. 'NO') then
    print *, 'PRECONDITIONER_FOR_SYSTEM_MATRIX    none'
  else if(answer .eq. 'DI') then
    print *, 'PRECONDITIONER_FOR_SYSTEM_MATRIX    diagonal'
  else if(answer .eq. 'IC') then
    print *, 'PRECONDITIONER_FOR_SYSTEM_MATRIX    incomplete_cholesky'
  endif

  print *, '# Tolerance for velocity solver: '
  call Tokenizer_Mod_Read_Line(CMN_FILE)
  print *, 'TOLERANCE_FOR_MOMENTUM_SOLVER    ', trim(line % tokens(1))
  print *, '# Tolerance for pressure solver: '
  call Tokenizer_Mod_Read_Line(CMN_FILE)
  print *, 'TOLERANCE_FOR_PRESSURE_SOLVER    ', trim(line % tokens(1))
  if(model .ne. 'LES' .and. model .ne. 'DNS') then
    print *, '# Tolerance for turbulence solver: '
    call Tokenizer_Mod_Read_Line(CMN_FILE)
    print *, 'TOLERANCE_FOR_TURBULENCE_SOLVER    ', trim(line % tokens(1))
  end if
  if(HOT) then
    print *, '# Tolerance for energy solver: '
    call Tokenizer_Mod_Read_Line(CMN_FILE)
    print *, 'TOLERANCE_FOR_ENERGY_SOLVER    ', trim(line % tokens(1))
  end if

  if(SIMPLE) then
    print *, '# Tolerance for SIMPLE: '
    call Tokenizer_Mod_Read_Line(CMN_FILE)
    print *, 'TOLERANCE_FOR_SIMPLE_ALGORITHM    ', trim(line % tokens(1))
  endif

  ! Time step
  print *, '# Time step: '
  call Tokenizer_Mod_Read_Line(CMN_FILE)
  print *, 'TIME_STEP    ', trim(line % tokens(1))

  ! Wall velocity
!@do m=1,grid % n_materials
    print *, '# Enter Pdrop (x, y, z) '
    call Tokenizer_Mod_Read_Line(CMN_FILE)
    print *, 'PRESSURE_DROPS    ', trim(line % tokens(1)), '  ',  &
                                   trim(line % tokens(2)), '  ',  &
                                   trim(line % tokens(3))
!@end do

  ! Mass fluxes
!@do m=1,grid % n_materials
    print *, '# Enter mass flow rates (x, y, z) '
    call Tokenizer_Mod_Read_Line(CMN_FILE)
    print *, 'MASS_FLOW_RATES    ', trim(line % tokens(1)), '  ',  &
                                    trim(line % tokens(2)), '  ',  &
                                    trim(line % tokens(3))
!@end do

  !--------------------------------------!
  !   Interpolate between diff. meshes   !
  !--------------------------------------!
  ! call Load_Ini(grid)

  !--------------------------------------------!
  !   Loading data from previous computation   !
  !--------------------------------------------!
  ! call Load_Restart_Ini(grid)

  !----------------------------!
  !   Close the command file   !
  !----------------------------!
  close(CMN_FILE)

  end program

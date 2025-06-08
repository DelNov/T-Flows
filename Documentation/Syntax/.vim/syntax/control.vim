" Vim syntax file
" Language:	control

au BufNewFile,BufRead control* setf control

" quit when a syntax file was already loaded
if exists("b:current_syntax")
   finish
endif

" ==============================================================================
" Comments:
syn region controlComment  start="\#"  end="\n"
syn region controlComment  start="\!"  end="\n"
syn region controlComment  start="%"   end="\n"

" List of keyords was obtained with command:
" cat ../Documentation/all_control_keywords | awk '{print "syn keyword controlKeyword   " $1}' | grep -v '|' | grep -v KEYWORD
" ==============================================================================
" Mandatory problem name and boundary conditions
syn keyword controlKeyword   PROBLEM_NAME
syn keyword controlKeyword   BOUNDARY_CONDITION
syn keyword controlKeyword   TYPE
syn keyword controlKeyword   VARIABLES
syn keyword controlKeyword   VALUES
syn keyword controlKeyword   FILE
syn keyword controlKeyword   INTERFACE_CONDITION
syn keyword controlKeyword   BOUNDARY_CONDITIONS
" Some other non-mandatory but often used kewords
syn keyword controlKeyword   INITIAL_CONDITION
syn keyword controlKeyword   SYNTHETIC_EDDIES
syn keyword controlKeyword   NUMBER_OF_EDDIES
syn keyword controlKeyword   MAX_EDDY_RADIUS
syn keyword controlKeyword   EDDY_INTENSITY
" Input-output
syn keyword controlKeyword   BACKUP_SAVE_INTERVAL
syn keyword controlKeyword   LOAD_BACKUP_NAME
syn keyword controlKeyword   NUMBER_OF_MONITORING_POINTS
syn keyword controlKeyword   MONITORING_POINT_001
syn keyword controlKeyword   MONITORING_POINT_002
syn keyword controlKeyword   MONITORING_POINT_003
syn keyword controlKeyword   MONITORING_POINT_004
syn keyword controlKeyword   MONITORING_POINT_005
syn keyword controlKeyword   MONITORING_POINT_006
syn keyword controlKeyword   MONITORING_POINT_007
syn keyword controlKeyword   MONITORING_POINT_008
syn keyword controlKeyword   MONITORING_POINT_009
syn keyword controlKeyword   MONITORING_POINT_010
syn keyword controlKeyword   MONITORING_POINT_011
syn keyword controlKeyword   MONITORING_POINT_012
syn keyword controlKeyword   MONITORING_POINT_013
syn keyword controlKeyword   MONITORING_POINT_014
syn keyword controlKeyword   MONITORING_POINT_015
syn keyword controlKeyword   MONITORING_POINT_016
syn keyword controlKeyword   RESULTS_SAVE_INTERVAL
syn keyword controlKeyword   SAVE_BACKUP_NAME
syn keyword controlKeyword   SAVE_INITIAL_CONDITION
syn keyword controlKeyword   SAVE_RESULTS_AT_BOUNDARIES
syn keyword controlKeyword   SAVE_RESULTS_UNITS  
syn keyword controlKeyword   SWARM_SAVE_INTERVAL
syn keyword controlKeyword   WALL_TIME_MAX_HOURS
syn keyword controlKeyword   PROFILER_INFO
" Numerics
syn keyword controlKeyword   ADVECTION_SCHEME_FOR_ENERGY
syn keyword controlKeyword   ADVECTION_SCHEME_FOR_MOMENTUM
syn keyword controlKeyword   ADVECTION_SCHEME_FOR_SCALARS
syn keyword controlKeyword   ADVECTION_SCHEME_FOR_TURBULENCE
syn keyword controlKeyword   ADVECTION_SCHEME_FOR_VOF
syn keyword controlKeyword   BLENDING_COEFFICIENT_FOR_ENERGY
syn keyword controlKeyword   BLENDING_COEFFICIENT_FOR_MOMENTUM
syn keyword controlKeyword   BLENDING_COEFFICIENT_FOR_SCALARS
syn keyword controlKeyword   BLENDING_COEFFICIENT_FOR_TURBULENCE
syn keyword controlKeyword   BLENDING_COEFFICIENT_FOR_VOF
syn keyword controlKeyword   BLEND_SYSTEM_MATRICES
syn keyword controlKeyword   CHOI_CORRECTION
syn keyword controlKeyword   GRADIENT_METHOD_FOR_ENERGY
syn keyword controlKeyword   GRADIENT_METHOD_FOR_MOMENTUM
syn keyword controlKeyword   GRADIENT_METHOD_FOR_PRESSURE
syn keyword controlKeyword   GRADIENT_METHOD_FOR_SCALARS
syn keyword controlKeyword   GRADIENT_METHOD_FOR_TURBULENCE
syn keyword controlKeyword   GRADIENT_METHOD_FOR_VOF
syn keyword controlKeyword   GRADIENT_METHOD_FOR_WALL_DISTANCE
syn keyword controlKeyword   GU_CORRECTION
syn keyword controlKeyword   LINEAR_SOLVERS
syn keyword controlKeyword   MAX_CORRECTION_CYCLES_BETA_VOF
syn keyword controlKeyword   MAX_COURANT_VOF
syn keyword controlKeyword   MAX_GAUSS_GRADIENTS_ITERATIONS
syn keyword controlKeyword   MAX_ITERATIONS_FOR_ENERGY_SOLVER
syn keyword controlKeyword   MAX_ITERATIONS_FOR_MOMENTUM_SOLVER
syn keyword controlKeyword   MAX_ITERATIONS_FOR_POTENTIAL_SOLVER
syn keyword controlKeyword   MAX_ITERATIONS_FOR_PRESSURE_SOLVER
syn keyword controlKeyword   MAX_ITERATIONS_FOR_SCALARS_SOLVER
syn keyword controlKeyword   MAX_ITERATIONS_FOR_TURBULENCE_SOLVER
syn keyword controlKeyword   MAX_ITERATIONS_FOR_VOF_SOLVER
syn keyword controlKeyword   MAX_ITERATIONS_FOR_WALL_DISTANCE_SOLVER
syn keyword controlKeyword   MAX_LEAST_SQUARES_GRADIENTS_ITERATIONS
syn keyword controlKeyword   MAX_SIMPLE_ITERATIONS
syn keyword controlKeyword   MAX_SMOOTH_CYCLES_CURVATURE_VOF
syn keyword controlKeyword   MAX_SUBSTEP_CYCLES_VOF
syn keyword controlKeyword   MAX_THREADS
syn keyword controlKeyword   MIN_SIMPLE_ITERATIONS
syn keyword controlKeyword   NORMALIZATION_FOR_ENERGY_SOLVER
syn keyword controlKeyword   NORMALIZATION_FOR_MOMENTUM_SOLVER
syn keyword controlKeyword   NORMALIZATION_FOR_PRESSURE_SOLVER
syn keyword controlKeyword   NORMALIZATION_FOR_SIMPLE_ALGORITHM
syn keyword controlKeyword   NORMALIZATION_FOR_SCALARS_SOLVER
syn keyword controlKeyword   NORMALIZATION_FOR_TURBULENCE_SOLVER
syn keyword controlKeyword   NUMBER_OF_PISO_CORRECTIONS
syn keyword controlKeyword   NUMBER_OF_TIME_STEPS
syn keyword controlKeyword   NUMBER_OF_SWARM_SUB_STEPS
syn keyword controlKeyword   PRECONDITIONER_FOR_SYSTEM_MATRIX
syn keyword controlKeyword   PRESSURE_MOMENTUM_COUPLING
syn keyword controlKeyword   REPORT_VOLUME_BALANCE
syn keyword controlKeyword   SIMPLE_UNDERRELAXATION_FOR_ENERGY
syn keyword controlKeyword   SIMPLE_UNDERRELAXATION_FOR_MOMENTUM
syn keyword controlKeyword   SIMPLE_UNDERRELAXATION_FOR_PRESSURE
syn keyword controlKeyword   SIMPLE_UNDERRELAXATION_FOR_SCALARS
syn keyword controlKeyword   SIMPLE_UNDERRELAXATION_FOR_TURBULENCE
syn keyword controlKeyword   SIMPLE_UNDERRELAXATION_FOR_VOF
syn keyword controlKeyword   SKEWNESS_CORRECTION_VOF
syn keyword controlKeyword   SOLVER_FOR_ENERGY
syn keyword controlKeyword   SOLVER_FOR_MOMENTUM
syn keyword controlKeyword   SOLVER_FOR_POTENTIAL
syn keyword controlKeyword   SOLVER_FOR_PRESSURE
syn keyword controlKeyword   SOLVER_FOR_SCALARS
syn keyword controlKeyword   SOLVER_FOR_TURBULENCE
syn keyword controlKeyword   SOLVER_FOR_VOF
syn keyword controlKeyword   SOLVER_FOR_WALL_DISTANCE
syn keyword controlKeyword   TIME_INTEGRATION_SCHEME
syn keyword controlKeyword   TIME_STEP
syn keyword controlKeyword   TOLERANCE_FOR_GAUSS_GRADIENTS
syn keyword controlKeyword   TOLERANCE_FOR_ENERGY_SOLVER
syn keyword controlKeyword   TOLERANCE_FOR_MOMENTUM_SOLVER
syn keyword controlKeyword   TOLERANCE_FOR_POTENTIAL_SOLVER
syn keyword controlKeyword   TOLERANCE_FOR_PRESSURE_SOLVER
syn keyword controlKeyword   TOLERANCE_FOR_SIMPLE_ALGORITHM
syn keyword controlKeyword   TOLERANCE_FOR_SCALARS_SOLVER
syn keyword controlKeyword   TOLERANCE_FOR_TURBULENCE_SOLVER
syn keyword controlKeyword   TOLERANCE_FOR_VOF_SOLVER
syn keyword controlKeyword   TOLERANCE_FOR_WALL_DISTANCE_SOLVER
" PETSc options
syn keyword controlKeyword   PETSC_OPTIONS
syn keyword controlKeyword   PETSC_OPTIONS_FOR_MOMENTUM
syn keyword controlKeyword   PETSC_OPTIONS_FOR_PRESSURE
syn keyword controlKeyword   PETSC_OPTIONS_FOR_WALL_DISTANCE
syn keyword controlKeyword   PETSC_OPTIONS_FOR_POTENTIAL
syn keyword controlKeyword   PETSC_OPTIONS_FOR_VOF
syn keyword controlKeyword   PETSC_OPTIONS_FOR_ENERGY
syn keyword controlKeyword   PETSC_OPTIONS_FOR_SCALARS
syn keyword controlKeyword   PETSC_OPTIONS_FOR_TURBULENCE
syn keyword controlKeyword   SOLVER
syn keyword controlKeyword   PREC
syn keyword controlKeyword   PREC_OPTS
syn keyword controlKeyword   TOLERANCE
" Physics
syn keyword controlKeyword   ANGULAR_VELOCITY_VECTOR
syn keyword controlKeyword   BULK_VELOCITIES
syn keyword controlKeyword   BUOYANCY
syn keyword controlKeyword   DYNAMIC_VISCOSITY
syn keyword controlKeyword   EXTRAPOLATE_TEMPERATURE_EXP
syn keyword controlKeyword   GRAVITATIONAL_VECTOR
syn keyword controlKeyword   HEAT_CAPACITY
syn keyword controlKeyword   HEAT_TRANSFER
syn keyword controlKeyword   HYBRID_LES_RANS_SWITCH
syn keyword controlKeyword   INTERFACE_TRACKING
syn keyword controlKeyword   LATENT_HEAT
syn keyword controlKeyword   LEE_MODEL_COEFFICIENTS
syn keyword controlKeyword   MASS_DENSITY
syn keyword controlKeyword   MASS_TRANSFER_MODEL
syn keyword controlKeyword   MAX_PARTICLES
syn keyword controlKeyword   NUMBER_OF_DOMAINS
syn keyword controlKeyword   NUMBER_OF_PHASES
syn keyword controlKeyword   NUMBER_OF_SCALARS
syn keyword controlKeyword   NUMBER_OF_SWARM_SUBSTEPS
syn keyword controlKeyword   PARTICLE_TRACKING
syn keyword controlKeyword   PHASE_DENSITIES
syn keyword controlKeyword   PHASE_VISCOSITIES
syn keyword controlKeyword   PHASE_CAPACITIES
syn keyword controlKeyword   PHASE_CONDUCTIVITIES
syn keyword controlKeyword   POINT_FOR_MONITORING_PLANES
syn keyword controlKeyword   POTENTIAL_INITIALIZATION
syn keyword controlKeyword   PRESSURE_DROPS
syn keyword controlKeyword   REFERENCE_DENSITY
syn keyword controlKeyword   REFERENCE_TEMPERATURE
syn keyword controlKeyword   ROUGHNESS_COEFFICIENT
syn keyword controlKeyword   ROUGH_WALLS
syn keyword controlKeyword   SATURATION_TEMPERATURE
syn keyword controlKeyword   SCALARS_DIFFUSIVITY
syn keyword controlKeyword   SMAGORINSKY_CONSTANT
syn keyword controlKeyword   STARTING_TIME_STEP_FOR_SWARM_COMPUTATION
syn keyword controlKeyword   STARTING_TIME_STEP_FOR_SWARM_STATISTICS
syn keyword controlKeyword   STARTING_TIME_STEP_FOR_TURB_STATISTICS
syn keyword controlKeyword   SURFACE_TENSION
syn keyword controlKeyword   SWARM_COEFFICIENT_OF_RESTITUTION
syn keyword controlKeyword   SWARM_DENSITY
syn keyword controlKeyword   SWARM_DIAMETER
syn keyword controlKeyword   SWARM_SUBGRID_SCALE_MODEL
syn keyword controlKeyword   THERMAL_CONDUCTIVITY
syn keyword controlKeyword   TRACK_FRONT
syn keyword controlKeyword   TRACK_SURFACE
syn keyword controlKeyword   TURBULENCE_MODEL
syn keyword controlKeyword   TURBULENCE_MODEL_VARIANT
syn keyword controlKeyword   TURBULENT_HEAT_FLUX_MODEL
syn keyword controlKeyword   TURBULENT_PRANDTL_NUMBER
syn keyword controlKeyword   TURBULENT_SCHMIDT_NUMBER
syn keyword controlKeyword   VOLUME_EXPANSION_COEFFICIENT
syn keyword controlKeyword   VOLUME_FLOW_RATES
" Porous regions
syn keyword controlKeyword   NUMBER_OF_POROUS_REGIONS
syn keyword controlKeyword   POROUS_REGION_001
syn keyword controlKeyword   POROUS_REGION_002
syn keyword controlKeyword   POROUS_REGION_003
syn keyword controlKeyword   POROUS_REGION_004
syn keyword controlKeyword   POROUS_REGION_005
syn keyword controlKeyword   POROUS_REGION_006
syn keyword controlKeyword   POROUS_REGION_007
syn keyword controlKeyword   POROUS_REGION_008
syn keyword controlKeyword   POROUS_REGION_009
syn keyword controlKeyword   POROUS_REGION_010
syn keyword controlKeyword   POROUS_REGION_011
syn keyword controlKeyword   POROUS_REGION_012
syn keyword controlKeyword   POROUS_REGION_013
syn keyword controlKeyword   POROUS_REGION_014
syn keyword controlKeyword   POROUS_REGION_015
syn keyword controlKeyword   POROUS_REGION_016
syn keyword controlKeyword   STL_FILE  FROM_GRID

" ==============================================================================
" Numbers (integer must be before the float, otherwise things get messed up)
syn match controlInteger   "-\=\<[0-9]*\>"
syn match controlFloat     "-\=\d\+\.\d*\([eE][-+]\=\d\+\)\=[fl]\=\>"

" Intrinsic (these are T-Flows variables)
syn keyword controlIntrinsic             u  v  w  t  q  p  kin  eps  zeta  f22  vis  uu  vv  ww  uv  uw  vw  vof  vof_c_ang
syn keyword controlIntrinsic             c_01  c_02  c_03  c_04  c_05  c_06  q_01  q_02  q_03  q_04  q_05  q_06
syn keyword controlIntrinsic             x  y  z  rx  ry  rz
syn keyword controlBoundaryConditon      wall  wall_flux  inflow  outflow  pressure  convective  symmetry
syn keyword controlLinearSolvers         native  petsc  cg  bicg  incomplete_cholesky  diagonal  none  asm  hypre  rs_amg
syn keyword controlLinearSolvers         log  log_view
syn keyword controlNumericalParameters   simple  piso  linear  parabolic  gauss_theorem  least_squares
syn keyword controlNumericalParameters   central  smart  luds  quick  smart  gamma  minmod  blended  superbee  avl_smart
syn keyword controlNumericalParameters   cicsam  upwind  stacs
syn keyword controlNumericalParameters   yes  no
syn keyword controlNumericalParameters   percentage  seconds
syn keyword controlPhysicalModels        thermal  density
syn keyword controlTurbulenceModels      k_eps_zeta_f  k_eps  les_tvm  les_wale  les_dynamic  les_smagorinsky  hybrid_les_rans  hybrid_les_prandtl
syn keyword controlTurbulenceModels      des_spalart  spalart_allmaras  dns  rsm_hanjalic_jakirlic  rsm_manceau_hanjalic  dns  none
syn keyword controlTurbulenceModels      high_re  low_re  ggdh  sgdh  afm
syn keyword controlMassTransferModels    temperature_gradients  lee

" ==============================================================================
" The default methods for highlighting. Can be overridden later.
hi def link controlKeyword               Type
hi def link controlComment               Comment
hi def link controlInteger               Number
hi def link controlFloat                 Number
hi def link controlIntrinsic             Identifier
hi def link controlBoundaryConditon      Keyword
hi def link controlLinearSolvers         Keyword
hi def link controlNumericalParameters   Keyword
hi def link controlPhysicalModels        Keyword
hi def link controlTurbulenceModels      Keyword
hi def link controlMassTransferModels    Keyword

let b:current_syntax = "control"


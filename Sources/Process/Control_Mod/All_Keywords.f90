    character(len=40), parameter :: ALL_KEYWORDS(78) = (/                      &
    !------------------!
    !   Input/Output   !
    !------------------!
    'BACKUP_SAVE_INTERVAL                    ',  &
    'LOAD_BACKUP_NAME                        ',  &
    'LOAD_INITIAL_SOLUTION_NAME              ',  &
    'NUMBER_OF_MONITORING_POINTS             ',  &
    'PROBLEM_NAME                            ',  &
    'RESULTS_SAVE_INTERVAL                   ',  &
    'SAVE_BACKUP_NAME                        ',  &
    'SAVE_INITIAL_SOLUTION_NAME              ',  &
    !--------------!
    !   Numerics   !
    !--------------!
    'ADVECTION_SCHEME_FOR_ENERGY             ',  &
    'ADVECTION_SCHEME_FOR_MOMENTUM           ',  &
    'ADVECTION_SCHEME_FOR_TURBULENCE         ',  &
    'BLENDING_COEFFICIENT_FOR_ENERGY         ',  &
    'BLENDING_COEFFICIENT_FOR_MOMENTUM       ',  &
    'BLENDING_COEFFICIENT_FOR_TURBULENCE     ',  &
    'MAX_ITERATIONS_FOR_ENERGY_SOLVER        ',  &
    'MAX_ITERATIONS_FOR_MOMENTUM_SOLVER      ',  &
    'MAX_ITERATIONS_FOR_PRESSURE_SOLVER      ',  &
    'MAX_ITERATIONS_FOR_TURBULENCE_SOLVER    ',  &
    'MAX_SIMPLE_ITERATIONS                   ',  &
    'MIN_SIMPLE_ITERATIONS                   ',  &
    'NORMALIZATION_FOR_ENERGY_SOLVER         ',  &
    'NORMALIZATION_FOR_MOMENTUM_SOLVER       ',  &
    'NORMALIZATION_FOR_PRESSURE_SOLVER       ',  &
    'NORMALIZATION_FOR_SIMPLE_ALGORITHM      ',  &
    'NORMALIZATION_FOR_TURBULENCE_SOLVER     ',  &
    'NUMBER_OF_TIME_STEPS                    ',  &
    'PRECONDITIONER_FOR_SYSTEM_MATRIX        ',  &
    'PRESSURE_MOMENTUM_COUPLING              ',  &
    'SIMPLE_UNDERRELAXATION_FOR_ENERGY       ',  &
    'SIMPLE_UNDERRELAXATION_FOR_MOMENTUM     ',  &
    'SIMPLE_UNDERRELAXATION_FOR_PRESSURE     ',  &
    'SIMPLE_UNDERRELAXATION_FOR_TURBULENCE   ',  &
    'SOLVER_FOR_ENERGY                       ',  &
    'SOLVER_FOR_MOMENTUM                     ',  &
    'SOLVER_FOR_PRESSURE                     ',  &
    'SOLVER_FOR_TURBULENCE                   ',  &
    'STARTING_TIME_STEP_FOR_STATISTICS       ',  &
    'TIME_INTEGRATION_FOR_ADVECTION          ',  &
    'TIME_INTEGRATION_FOR_CROSS_DIFFUSION    ',  &
    'TIME_INTEGRATION_FOR_DIFFUSION          ',  &
    'TIME_INTEGRATION_FOR_INERTIA            ',  &
    'TIME_STEP                               ',  &
    'TOLERANCE_FOR_ENERGY_SOLVER             ',  &
    'TOLERANCE_FOR_MOMENTUM_SOLVER           ',  &
    'TOLERANCE_FOR_PRESSURE_SOLVER           ',  &
    'TOLERANCE_FOR_SIMPLE_ALGORITHM          ',  &
    'TOLERANCE_FOR_TURBULENCE_SOLVER         ',  &
    !-------------!
    !   Physics   !
    !-------------!
    'ANGULAR_VELOCITY_VECTOR                 ',  &
    'BUOYANCY                                ',  &
    'DYNAMIC_VISCOSITY                       ',  &
    'GRAVITATIONAL_VECTOR                    ',  &
    'HEAT_CAPACITY                           ',  &
    'HEAT_TRANSFER                           ',  &
    'MASS_DENSITY                            ',  &
    'MASS_FLOW_RATES                         ',  &
    'NUMBER_OF_PHASES                        ',  &
    'NUMBER_OF_SPECIES                       ',  &
    'POINT_FOR_MONITORING_PLANES             ',  &
    'PRESSURE_DROPS                          ',  &
    'REFERENCE_TEMPERATURE                   ',  &
    'ROUGHNESS_COEFFICIENT                   ',  &
    'ROUGH_WALLS                             ',  &
    'THERMAL_CONDUCTIVITY                    ',  &
    'TURBULENCE_MODEL                        ',  &
    'TURBULENCE_MODEL_VARIANT                ',  &
    'TURBULENCE_STATISTICS                   ',  &
    'TURBULENCE_WALL_TREATMENT               ',  &
    'TURBULENT_PRANDTL_NUMBER                ',  &
    'TURBULENT_SCHMIDT_NUMBER                ',  &
    !-------------------------------------------!
    !   User (should Schmidt number be here?)   !
    !-------------------------------------------!
    'ADVECTION_SCHEME_FOR_USER_SCALARS       ',  &
    'BLENDING_COEFFICIENT_FOR_USER_SCALARS   ',  &
    'MAX_ITERATIONS_FOR_USER_SCALARS_SOLVER  ',  &
    'NORMALIZATION_FOR_USER_SCALARS_SOLVER   ',  &
    'NUMBER_OF_USER_ARRAYS                   ',  &
    'NUMBER_OF_USER_SCALARS                  ',  &
    'SIMPLE_UNDERRELAXATION_FOR_USER_SCALARS ',  &
    'SOLVER_FOR_USER_SCALARS                 ',  &
    'TOLERANCE_FOR_USER_SCALARS_SOLVER       '   &
    /)

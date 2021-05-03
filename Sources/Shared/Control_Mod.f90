!==============================================================================!
  module Control_Mod
!----------------------------------[Modules]-----------------------------------!
  use Math_Mod
  use File_Mod
  use Comm_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Control file unit
  integer :: control_file_unit          ! current control file unit
  integer :: root_control_file_unit     ! root control file unit
  integer :: dom_control_file_unit(MD)  ! domain control file units

  ! Number of presumably misstyped keywords followed by their list
  integer       :: n_similar =  0
  character(SL) :: similar(128)

  contains

  !--------------------!
  !   Input / Output   !
  !--------------------!

    ! Basic functionality
    include 'Control_Mod/Basic_Functions/Open_Domain_File.f90'
    include 'Control_Mod/Basic_Functions/Open_Root_File.f90'
    include 'Control_Mod/Basic_Functions/Position_At_One_Key.f90'
    include 'Control_Mod/Basic_Functions/Position_At_Two_Keys.f90'
    include 'Control_Mod/Basic_Functions/Position_At_Three_Keys.f90'
    include 'Control_Mod/Basic_Functions/Read_Char_Item.f90'
    include 'Control_Mod/Basic_Functions/Read_Char_Item_On.f90'
    include 'Control_Mod/Basic_Functions/Read_Strings_On.f90'
    include 'Control_Mod/Basic_Functions/Read_Int_Item.f90'
    include 'Control_Mod/Basic_Functions/Read_Int_Item_On.f90'
    include 'Control_Mod/Basic_Functions/Read_Real_Item.f90'
    include 'Control_Mod/Basic_Functions/Read_Real_Item_On.f90'
    include 'Control_Mod/Basic_Functions/Read_Real_Array.f90'
    include 'Control_Mod/Basic_Functions/Read_Real_Array_On.f90'
    include 'Control_Mod/Basic_Functions/Similar_Warning.f90'
    include 'Control_Mod/Basic_Functions/Switch_To_Domain.f90'
    include 'Control_Mod/Basic_Functions/Switch_To_Root.f90'
!   include 'Control_Mod/Basic_Functions/Write_File.f90'

    ! Load
    include 'Control_Mod/Input_Output/Problem_Name.f90'
    include 'Control_Mod/Input_Output/Load_Backup_Name.f90'
    include 'Control_Mod/Input_Output/Save_Backup_Name.f90'
    include 'Control_Mod/Input_Output/Save_Initial_Condition.f90'
    include 'Control_Mod/Input_Output/Load_Initial_Solution_Name.f90'
    include 'Control_Mod/Input_Output/Save_Initial_Solution_Name.f90'
    include 'Control_Mod/Input_Output/Wall_Time_Max_Hours.f90'

    ! Save
    include 'Control_Mod/Input_Output/Backup_Save_Interval.f90'
    include 'Control_Mod/Input_Output/Results_Save_Interval.f90'
    include 'Control_Mod/Input_Output/Swarm_Save_Interval.f90'

  !-------------!
  !   Physics   !
  !-------------!

    ! Number of domains you are solving
    include 'Control_Mod/Physics/Number_Of_Domains.f90'

    ! Physical properties
    include 'Control_Mod/Physics/Dynamic_Viscosity.f90'
    include 'Control_Mod/Physics/Heat_Capacity.f90'
    include 'Control_Mod/Physics/Latent_Heat.f90'
    include 'Control_Mod/Physics/Mass_Density.f90'
    include 'Control_Mod/Physics/Thermal_Conductivity.f90'
    include 'Control_Mod/Physics/Scalars_Diffusivity.f90'

    ! Multiphase flow
    include 'Control_Mod/Physics/Multiphase_Model.f90'
    include 'Control_Mod/Physics/Mass_Transfer.f90'

    ! Heat transfer
    include 'Control_Mod/Physics/Heat_Transfer.f90'
    include 'Control_Mod/Physics/Buoyancy.f90'
    include 'Control_Mod/Physics/Reference_Temperature.f90'
    include 'Control_Mod/Physics/Saturation_Temperature.f90'
    include 'Control_Mod/Physics/Volume_Expansion_Coefficient.f90'

    ! Turbulence
    include 'Control_Mod/Physics/Roughness_Coefficient.f90'
    include 'Control_Mod/Physics/Rough_Walls.f90'
    include 'Control_Mod/Physics/Smagorinsky_Constant.f90'
    include 'Control_Mod/Physics/Turbulence_Model.f90'
    include 'Control_Mod/Physics/Turbulence_Model_Variant.f90'
    include 'Control_Mod/Physics/Turbulent_Prandtl_Number.f90'
    include 'Control_Mod/Physics/Turbulent_Schmidt_Number.f90'
    include 'Control_Mod/Physics/Turbulent_Heat_Flux_Model.f90'
    include 'Control_Mod/Physics/Hybrid_Les_Rans_Switch.f90'

    ! Other environmental conditions
    include 'Control_Mod/Physics/Angular_Velocity_Vector.f90'
    include 'Control_Mod/Physics/Gravitational_Vector.f90'
    include 'Control_Mod/Physics/Mass_Flow_Rates.f90'
    include 'Control_Mod/Physics/Pressure_Drops.f90'
    include 'Control_Mod/Physics/Point_For_Monitoring_Planes.f90'
    include 'Control_Mod/Physics/Potential_Initialization.f90'

    ! Multiphase
    include 'Control_Mod/Physics/Number_Of_Phases.f90'
    include 'Control_Mod/Physics/Number_Of_Particles.f90'
    include 'Control_Mod/Physics/Number_Of_Swarm_Sub_Steps.f90'
    include 'Control_Mod/Physics/Swarm_Coefficient_Of_Restitution.f90'
    include 'Control_Mod/Physics/Swarm_Density.f90'
    include 'Control_Mod/Physics/Swarm_Diameter.f90'
    include 'Control_Mod/Physics/Swarm_Subgrid_Scale_Model.f90'
    include 'Control_Mod/Physics/Phase_Capacities.f90'
    include 'Control_Mod/Physics/Phase_Conductivities.f90'
    include 'Control_Mod/Physics/Phase_Densities.f90'
    include 'Control_Mod/Physics/Phase_Viscosities.f90'
    include 'Control_Mod/Physics/Reference_Density.f90'
    include 'Control_Mod/Physics/Surface_Tension.f90'
    include 'Control_Mod/Physics/Track_Front.f90'

    ! Scalars (like species, for example)
    include 'Control_Mod/Physics/Number_Of_Scalars.f90'

    ! Statistics
    include 'Control_Mod/Physics/Starting_Time_Step_For_Turb_Statistics.f90'
    include 'Control_Mod/Physics/Starting_Time_Step_For_Swarm_Computation.f90'
    include 'Control_Mod/Physics/Starting_Time_Step_For_Swarm_Statistics.f90'

  !--------------!
  !   Numerics   !
  !--------------!

    ! Time Stepping
    include 'Control_Mod/Numerics/Time_Step.f90'
    include 'Control_Mod/Numerics/Number_Of_Time_Steps.f90'

    ! Discretization
    include 'Control_Mod/Numerics/Advection_Scheme_For_Energy.f90'
    include 'Control_Mod/Numerics/Advection_Scheme_For_Momentum.f90'
    include 'Control_Mod/Numerics/Advection_Scheme_For_Scalars.f90'
    include 'Control_Mod/Numerics/Advection_Scheme_For_Turbulence.f90'
    include 'Control_Mod/Numerics/Advection_Scheme_For_Multiphase.f90'
    include 'Control_Mod/Numerics/Blending_Coefficient_For_Energy.f90'
    include 'Control_Mod/Numerics/Blending_Coefficient_For_Momentum.f90'
    include 'Control_Mod/Numerics/Blending_Coefficient_For_Scalars.f90'
    include 'Control_Mod/Numerics/Blending_Coefficient_For_Turbulence.f90'
    include 'Control_Mod/Numerics/Blending_Coefficient_For_Multiphase.f90'
    include 'Control_Mod/Numerics/Choi_Correction.f90'
    include 'Control_Mod/Numerics/Gradient_Method_For_Energy.f90'
    include 'Control_Mod/Numerics/Gradient_Method_For_Momentum.f90'
    include 'Control_Mod/Numerics/Gradient_Method_For_Pressure.f90'
    include 'Control_Mod/Numerics/Gradient_Method_For_Scalars.f90'
    include 'Control_Mod/Numerics/Gradient_Method_For_Turbulence.f90'
    include 'Control_Mod/Numerics/Gradient_Method_For_Multiphase.f90'
    include 'Control_Mod/Numerics/Max_Iterations_For_Gauss_Gradients.f90'
    include 'Control_Mod/Numerics/Max_Simple_Iterations.f90'
    include 'Control_Mod/Numerics/Min_Simple_Iterations.f90'
    include 'Control_Mod/Numerics/Number_Of_Piso_Corrections.f90'
    include 'Control_Mod/Numerics/Pressure_Momentum_Coupling.f90'
    include 'Control_Mod/Numerics/Simple_Underrelaxation_For_Momentum.f90'
    include 'Control_Mod/Numerics/Simple_Underrelaxation_For_Pressure.f90'
    include 'Control_Mod/Numerics/Simple_Underrelaxation_For_Energy.f90'
    include 'Control_Mod/Numerics/Simple_Underrelaxation_For_Scalars.f90'
    include 'Control_Mod/Numerics/Simple_Underrelaxation_For_Turbulence.f90'
    include 'Control_Mod/Numerics/Simple_Underrelaxation_For_Multiphase.f90'
    include 'Control_Mod/Numerics/Time_Integration_Scheme.f90'
    include 'Control_Mod/Numerics/Tolerance_For_Gauss_Gradients.f90'

    ! Linear solvers
    include 'Control_Mod/Numerics/Preconditioner_For_System_Matrix.f90'
    include 'Control_Mod/Numerics/Solver_For_Energy.f90'
    include 'Control_Mod/Numerics/Solver_For_Momentum.f90'
    include 'Control_Mod/Numerics/Solver_For_Pressure.f90'
    include 'Control_Mod/Numerics/Solver_For_Scalars.f90'
    include 'Control_Mod/Numerics/Solver_For_Turbulence.f90'
    include 'Control_Mod/Numerics/Solver_For_Multiphase.f90'
    include 'Control_Mod/Numerics/Max_Iterations_For_Energy_Solver.f90'
    include 'Control_Mod/Numerics/Max_Iterations_For_Momentum_Solver.f90'
    include 'Control_Mod/Numerics/Max_Iterations_For_Presssure_Solver.f90'
    include 'Control_Mod/Numerics/Max_Iterations_For_Scalars_Solver.f90'
    include 'Control_Mod/Numerics/Max_Iterations_For_Turbulence_Solver.f90'
    include 'Control_Mod/Numerics/Max_Iterations_For_Multiphase_Solver.f90'
    include 'Control_Mod/Numerics/Normalization_For_Momentum_Solver.f90'
    include 'Control_Mod/Numerics/Normalization_For_Pressure_Solver.f90'
    include 'Control_Mod/Numerics/Normalization_For_Energy_Solver.f90'
    include 'Control_Mod/Numerics/Normalization_For_Scalars_Solver.f90'
    include 'Control_Mod/Numerics/Normalization_For_Simple_Algorithm.f90'
    include 'Control_Mod/Numerics/Normalization_For_Turbulence_Solver.f90'
    include 'Control_Mod/Numerics/Tolerance_For_Momentum_Solver.f90'
    include 'Control_Mod/Numerics/Tolerance_For_Pressure_Solver.f90'
    include 'Control_Mod/Numerics/Tolerance_For_Energy_Solver.f90'
    include 'Control_Mod/Numerics/Tolerance_For_Scalars_Solver.f90'
    include 'Control_Mod/Numerics/Tolerance_For_Turbulence_Solver.f90'
    include 'Control_Mod/Numerics/Tolerance_For_Multiphase_Solver.f90'
    include 'Control_Mod/Numerics/Tolerance_For_Simple_Algorithm.f90'
    include 'Control_Mod/Numerics/V_Cycle_Max_Grid_Levels.f90'
    include 'Control_Mod/Numerics/V_Cycle_Number_Of_Smoothing_Iterations.f90'
    include 'Control_Mod/Numerics/V_Cycle_Residual_Ratio.f90'

    ! Numerical Parameters VOF (CICSAM)
    include 'Control_Mod/Numerics/Max_Correction_Cycles_Beta_Vof.f90'
    include 'Control_Mod/Numerics/Max_Smoothing_Cycles_Curvature_Vof.f90'
    include 'Control_Mod/Numerics/Max_Smoothing_Cycles_Normal_Vof.f90'
    include 'Control_Mod/Numerics/Max_Courant_Vof.f90'
    include 'Control_Mod/Numerics/Max_Substep_Cycles_Vof.f90'
    include 'Control_Mod/Numerics/Skewness_Correction_Vof.f90'
  end module

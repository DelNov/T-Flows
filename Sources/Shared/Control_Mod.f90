#include "../Shared/Assert.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Control_Mod
!------------------------------------------------------------------------------!
!>  The Control_Mod module is a central component of controling a simulation
!>  in the T-Flows code.  It reads and processes the control file, which
!>  contains parameters and settings for a CFD simulations.
!----------------------------------[Modules]-----------------------------------!
  use Math_Mod
  use File_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !   Control type   !
  !------------------!
  !> Encapsulates various aspects of the control file, including file units
  !> (root control file and control files for different domains) and procedures
  !> for reading and manipulating control file data.  Procedures are classified
  !> in five big groups: basic functionality, input-output, native solvers,
  !> numerical models and physical models.  In most cases, each procedure in
  !> the Control_Type is reading one keyword in the control file (exceptions
  !> are the boundary and initial conditions).
  type Control_Type

    ! Control file unit
    integer, private :: file_unit          ! current control file unit
    integer, private :: root_file_unit     ! root control file unit
    integer, private :: dom_file_unit(MD)  ! domain control file units

    contains

      ! Basic functionality (manouvering through control file)
      procedure :: Open_Domain_File
      procedure :: Open_Root_File
      procedure :: Position_At_One_Key
      procedure :: Position_At_Two_Keys
      procedure :: Position_At_Three_Keys
      procedure :: Read_Char_Item
      procedure :: Read_Char_Item_On
      procedure :: Read_Char_Real_Vector
      procedure :: Read_Int_Item
      procedure :: Read_Int_Item_On
      procedure :: Read_Keyless_Strings_On
      procedure :: Read_Real_Item
      procedure :: Read_Real_Item_On
      procedure :: Read_Real_Vector
      procedure :: Read_Real_Vector_On
      procedure :: Read_Strings_On
      procedure :: Switch_To_Domain
      procedure :: Switch_To_Root

      ! Input/output
      procedure :: Read_Problem_Name
      procedure :: Load_Backup_Name
      procedure :: Save_Backup_Name
      procedure :: Save_Initial_Condition
      procedure :: Save_Results_At_Boundaries
      procedure :: Save_Results_Units
      procedure :: Wall_Time_Max_Hours
      procedure :: Backup_Save_Interval
      procedure :: Results_Save_Interval
      procedure :: Swarm_Save_Interval

      ! Unclassified as of yet
      procedure :: Linear_Solvers
      procedure :: Max_Threads

      ! Native solvers
      procedure :: Max_Iterations_For_Energy_Solver
      procedure :: Max_Iterations_For_Momentum_Solver
      procedure :: Max_Iterations_For_Potential_Solver
      procedure :: Max_Iterations_For_Pressure_Solver
      procedure :: Max_Iterations_For_Scalars_Solver
      procedure :: Max_Iterations_For_Turbulence_Solver
      procedure :: Max_Iterations_For_Vof_Solver
      procedure :: Max_Iterations_For_Wall_Distance_Solver
      procedure :: Preconditioner_For_System_Matrix
      procedure :: Solver_For_Energy
      procedure :: Solver_For_Momentum
      procedure :: Solver_For_Potential
      procedure :: Solver_For_Pressure
      procedure :: Solver_For_Scalars
      procedure :: Solver_For_Turbulence
      procedure :: Solver_For_Vof
      procedure :: Solver_For_Wall_Distance
      procedure :: Normalization_For_Momentum_Solver
      procedure :: Normalization_For_Pressure_Solver
      procedure :: Normalization_For_Energy_Solver
      procedure :: Normalization_For_Scalars_Solver
      procedure :: Normalization_For_Turbulence_Solver
      procedure :: Tolerance_For_Momentum_Solver
      procedure :: Tolerance_For_Potential_Solver
      procedure :: Tolerance_For_Pressure_Solver
      procedure :: Tolerance_For_Energy_Solver
      procedure :: Tolerance_For_Scalars_Solver
      procedure :: Tolerance_For_Turbulence_Solver
      procedure :: Tolerance_For_Vof_Solver
      procedure :: Tolerance_For_Wall_Distance_Solver

      ! Numerics
      procedure :: Time_Step
      procedure :: Number_Of_Time_Steps
      procedure :: Advection_Scheme_For_Energy
      procedure :: Advection_Scheme_For_Momentum
      procedure :: Advection_Scheme_For_Scalars
      procedure :: Advection_Scheme_For_Turbulence
      procedure :: Advection_Scheme_For_Vof
      procedure :: Blending_Coefficients_For_Energy
      procedure :: Blending_Coefficients_For_Momentum
      procedure :: Blending_Coefficients_For_Scalars
      procedure :: Blending_Coefficients_For_Turbulence
      procedure :: Blending_Coefficients_For_Vof
      procedure :: Blend_System_Matrices
      procedure :: Choi_Correction
      procedure :: Gradient_Method_For_Energy
      procedure :: Gradient_Method_For_Momentum
      procedure :: Gradient_Method_For_Pressure
      procedure :: Gradient_Method_For_Scalars
      procedure :: Gradient_Method_For_Turbulence
      procedure :: Gradient_Method_For_Vof
      procedure :: Gradient_Method_For_Wall_Distance
      procedure :: Gu_Correction
      procedure :: Max_Gauss_Gradients_Iterations
      procedure :: Max_Least_Squares_Gradients_Iterations
      procedure :: Max_Simple_Iterations
      procedure :: Min_Simple_Iterations
      procedure :: Number_Of_Piso_Corrections
      procedure :: Pressure_Momentum_Coupling
      procedure :: Report_Volume_Balance
      procedure :: Simple_Underrelaxation_For_Momentum
      procedure :: Simple_Underrelaxation_For_Pressure
      procedure :: Simple_Underrelaxation_For_Energy
      procedure :: Simple_Underrelaxation_For_Scalars
      procedure :: Simple_Underrelaxation_For_Turbulence
      procedure :: Simple_Underrelaxation_For_Vof
      procedure :: Time_Integration_Scheme
      procedure :: Tolerance_For_Gauss_Gradients
      procedure :: Tolerance_For_Simple_Algorithm
      procedure :: Max_Correction_Cycles_Beta_Vof
      procedure :: Max_Smooth_Cycles_Curvature_Vof
      procedure :: Max_Courant_Vof
      procedure :: Max_Substep_Cycles_Vof
      procedure :: Skewness_Correction_Vof

      ! Physics
      procedure :: Bulk_Velocities
      procedure :: Number_Of_Domains
      procedure :: Dynamic_Viscosity
      procedure :: Heat_Capacity
      procedure :: Latent_Heat
      procedure :: Mass_Density
      procedure :: Thermal_Conductivity
      procedure :: Scalars_Diffusivity
      procedure :: Scalars_Diffusivities
      procedure :: Interface_Tracking
      procedure :: Particle_Tracking
      procedure :: Mass_Transfer_Model
      procedure :: Lee_Model_Coefficients
      procedure :: Heat_Transfer
      procedure :: Buoyancy
      procedure :: Reference_Temperature
      procedure :: Saturation_Temperature
      procedure :: Volume_Expansion_Coefficient
      procedure :: Hybrid_Les_Rans_Switch
      procedure :: Roughness_Coefficient
      procedure :: Rough_Walls
      procedure :: Monin_Obukov
      procedure :: Smagorinsky_Constant
      procedure :: Wale_Constant
      procedure :: Turbulence_Model
      procedure :: Turbulent_Prandtl_Number
      procedure :: Turbulent_Schmidt_Number
      procedure :: Turbulent_Heat_Flux_Model
      procedure :: Turbulent_Scalar_Flux_Model
      procedure :: Gravitational_Vector
      procedure :: Volume_Flow_Rates
      procedure :: Pressure_Drops
      procedure :: Point_For_Monitoring_Planes
      procedure :: Potential_Initialization
      procedure :: Max_Particles
      procedure :: Number_Of_Phases
      procedure :: Number_Of_Swarm_Sub_Steps
      procedure :: Swarm_Coefficient_Of_Restitution
      procedure :: Swarm_Density
      procedure :: Swarm_Diameter
      procedure :: Swarm_Subgrid_Scale_Model
      procedure :: Phase_Capacities
      procedure :: Phase_Conductivities
      procedure :: Phase_Densities
      procedure :: Phase_Viscosities
      procedure :: Reference_Density
      procedure :: Surface_Tension
      procedure :: Track_Front
      procedure :: Track_Surface
      procedure :: Number_Of_Scalars
      procedure :: Starting_Time_Step_For_Turb_Statistics
      procedure :: Starting_Time_Step_For_Swarm_Computation
      procedure :: Starting_Time_Step_For_Swarm_Statistics
      procedure :: Extrapolate_Temperature_Exp

  end type

  type(Control_Type) :: Control  !! Singleton object of the class Control_Type
                                 !! to make the access to Control_Type
                                 !! procedures easier.
  contains

    !-------------------------!
    !   Basic functionality   !  (manouvering through control file)
    !-------------------------!

#   include "Control_Mod/Basic_Functions/Open_Domain_File.f90"
#   include "Control_Mod/Basic_Functions/Open_Root_File.f90"
#   include "Control_Mod/Basic_Functions/Position_At_One_Key.f90"
#   include "Control_Mod/Basic_Functions/Position_At_Two_Keys.f90"
#   include "Control_Mod/Basic_Functions/Position_At_Three_Keys.f90"
#   include "Control_Mod/Basic_Functions/Read_Char_Item.f90"
#   include "Control_Mod/Basic_Functions/Read_Char_Item_On.f90"
#   include "Control_Mod/Basic_Functions/Read_Char_Real_Vector.f90"
#   include "Control_Mod/Basic_Functions/Read_Int_Item.f90"
#   include "Control_Mod/Basic_Functions/Read_Int_Item_On.f90"
#   include "Control_Mod/Basic_Functions/Read_Keyless_Strings_On.f90"
#   include "Control_Mod/Basic_Functions/Read_Real_Item.f90"
#   include "Control_Mod/Basic_Functions/Read_Real_Item_On.f90"
#   include "Control_Mod/Basic_Functions/Read_Real_Vector.f90"
#   include "Control_Mod/Basic_Functions/Read_Real_Vector_On.f90"
#   include "Control_Mod/Basic_Functions/Read_Strings_On.f90"
#   include "Control_Mod/Basic_Functions/Switch_To_Domain.f90"
#   include "Control_Mod/Basic_Functions/Switch_To_Root.f90"

    !--------------------!
    !   Input / Output   !
    !--------------------!

    ! Load
#   include "Control_Mod/Input_Output/Problem_Name.f90"
#   include "Control_Mod/Input_Output/Load_Backup_Name.f90"
#   include "Control_Mod/Input_Output/Save_Backup_Name.f90"
#   include "Control_Mod/Input_Output/Save_Initial_Condition.f90"
#   include "Control_Mod/Input_Output/Save_Results_At_Boundaries.f90"
#   include "Control_Mod/Input_Output/Save_Results_Units.f90"
#   include "Control_Mod/Input_Output/Wall_Time_Max_Hours.f90"

    ! Save
#   include "Control_Mod/Input_Output/Backup_Save_Interval.f90"
#   include "Control_Mod/Input_Output/Results_Save_Interval.f90"
#   include "Control_Mod/Input_Output/Swarm_Save_Interval.f90"

    !--------------------!
    !   Linear solvers   !  (native (from T-Flows) or PETSc)
    !--------------------!

#   include "Control_Mod/Linear_Solvers.f90"

    !-------------!
    !   Threads   !
    !-------------!

#   include "Control_Mod/Max_Threads.f90"

    !--------------------!
    !   Native solvers   !
    !--------------------!

    ! Linear solvers
#   include "Control_Mod/Native/Max_Iterations_For_Energy_Solver.f90"
#   include "Control_Mod/Native/Max_Iterations_For_Momentum_Solver.f90"
#   include "Control_Mod/Native/Max_Iterations_For_Potential_Solver.f90"
#   include "Control_Mod/Native/Max_Iterations_For_Pressure_Solver.f90"
#   include "Control_Mod/Native/Max_Iterations_For_Scalars_Solver.f90"
#   include "Control_Mod/Native/Max_Iterations_For_Turbulence_Solver.f90"
#   include "Control_Mod/Native/Max_Iterations_For_Vof_Solver.f90"
#   include "Control_Mod/Native/Max_Iterations_For_Wall_Distance_Solver.f90"
#   include "Control_Mod/Native/Preconditioner_For_System_Matrix.f90"
#   include "Control_Mod/Native/Solver_For_Energy.f90"
#   include "Control_Mod/Native/Solver_For_Momentum.f90"
#   include "Control_Mod/Native/Solver_For_Potential.f90"
#   include "Control_Mod/Native/Solver_For_Pressure.f90"
#   include "Control_Mod/Native/Solver_For_Scalars.f90"
#   include "Control_Mod/Native/Solver_For_Turbulence.f90"
#   include "Control_Mod/Native/Solver_For_Vof.f90"
#   include "Control_Mod/Native/Solver_For_Wall_Distance.f90"
#   include "Control_Mod/Native/Normalization_For_Momentum_Solver.f90"
#   include "Control_Mod/Native/Normalization_For_Pressure_Solver.f90"
#   include "Control_Mod/Native/Normalization_For_Energy_Solver.f90"
#   include "Control_Mod/Native/Normalization_For_Scalars_Solver.f90"
#   include "Control_Mod/Native/Normalization_For_Turbulence_Solver.f90"
#   include "Control_Mod/Native/Tolerance_For_Momentum_Solver.f90"
#   include "Control_Mod/Native/Tolerance_For_Potential_Solver.f90"
#   include "Control_Mod/Native/Tolerance_For_Pressure_Solver.f90"
#   include "Control_Mod/Native/Tolerance_For_Energy_Solver.f90"
#   include "Control_Mod/Native/Tolerance_For_Scalars_Solver.f90"
#   include "Control_Mod/Native/Tolerance_For_Turbulence_Solver.f90"
#   include "Control_Mod/Native/Tolerance_For_Vof_Solver.f90"
#   include "Control_Mod/Native/Tolerance_For_Wall_Distance_Solver.f90"

    !--------------!
    !   Numerics   !
    !--------------!

    ! Time Stepping
#   include "Control_Mod/Numerics/Time_Step.f90"
#   include "Control_Mod/Numerics/Number_Of_Time_Steps.f90"

    ! Discretization
#   include "Control_Mod/Numerics/Advection_Scheme_For_Energy.f90"
#   include "Control_Mod/Numerics/Advection_Scheme_For_Momentum.f90"
#   include "Control_Mod/Numerics/Advection_Scheme_For_Scalars.f90"
#   include "Control_Mod/Numerics/Advection_Scheme_For_Turbulence.f90"
#   include "Control_Mod/Numerics/Advection_Scheme_For_Vof.f90"
#   include "Control_Mod/Numerics/Blending_Coefficients_For_Energy.f90"
#   include "Control_Mod/Numerics/Blending_Coefficients_For_Momentum.f90"
#   include "Control_Mod/Numerics/Blending_Coefficients_For_Scalars.f90"
#   include "Control_Mod/Numerics/Blending_Coefficients_For_Turbulence.f90"
#   include "Control_Mod/Numerics/Blending_Coefficients_For_Vof.f90"
#   include "Control_Mod/Numerics/Blend_System_Matrices.f90"
#   include "Control_Mod/Numerics/Choi_Correction.f90"
#   include "Control_Mod/Numerics/Gradient_Method_For_Energy.f90"
#   include "Control_Mod/Numerics/Gradient_Method_For_Momentum.f90"
#   include "Control_Mod/Numerics/Gradient_Method_For_Pressure.f90"
#   include "Control_Mod/Numerics/Gradient_Method_For_Scalars.f90"
#   include "Control_Mod/Numerics/Gradient_Method_For_Turbulence.f90"
#   include "Control_Mod/Numerics/Gradient_Method_For_Vof.f90"
#   include "Control_Mod/Numerics/Gradient_Method_For_Wall_Distance.f90"
#   include "Control_Mod/Numerics/Gu_Correction.f90"
#   include "Control_Mod/Numerics/Max_Gauss_Gradients_Iterations.f90"
#   include "Control_Mod/Numerics/Max_Least_Squares_Gradients_Iterations.f90"
#   include "Control_Mod/Numerics/Max_Simple_Iterations.f90"
#   include "Control_Mod/Numerics/Min_Simple_Iterations.f90"
#   include "Control_Mod/Numerics/Number_Of_Piso_Corrections.f90"
#   include "Control_Mod/Numerics/Pressure_Momentum_Coupling.f90"
#   include "Control_Mod/Numerics/Report_Volume_Balance.f90"
#   include "Control_Mod/Numerics/Simple_Underrelaxation_For_Momentum.f90"
#   include "Control_Mod/Numerics/Simple_Underrelaxation_For_Pressure.f90"
#   include "Control_Mod/Numerics/Simple_Underrelaxation_For_Energy.f90"
#   include "Control_Mod/Numerics/Simple_Underrelaxation_For_Scalars.f90"
#   include "Control_Mod/Numerics/Simple_Underrelaxation_For_Turbulence.f90"
#   include "Control_Mod/Numerics/Simple_Underrelaxation_For_Vof.f90"
#   include "Control_Mod/Numerics/Time_Integration_Scheme.f90"
#   include "Control_Mod/Numerics/Tolerance_For_Gauss_Gradients.f90"
#   include "Control_Mod/Numerics/Tolerance_For_Simple_Algorithm.f90"

    ! Numerical Parameters VOF (CICSAM)
#   include "Control_Mod/Numerics/Max_Correction_Cycles_Beta_Vof.f90"
#   include "Control_Mod/Numerics/Max_Smooth_Cycles_Curvature_Vof.f90"
#   include "Control_Mod/Numerics/Max_Courant_Vof.f90"
#   include "Control_Mod/Numerics/Max_Substep_Cycles_Vof.f90"
#   include "Control_Mod/Numerics/Skewness_Correction_Vof.f90"

    !-------------!
    !   Physics   !
    !-------------!

    ! Number of domains you are solving
#   include "Control_Mod/Physics/Number_Of_Domains.f90"

    ! Physical properties
#   include "Control_Mod/Physics/Dynamic_Viscosity.f90"
#   include "Control_Mod/Physics/Heat_Capacity.f90"
#   include "Control_Mod/Physics/Latent_Heat.f90"
#   include "Control_Mod/Physics/Mass_Density.f90"
#   include "Control_Mod/Physics/Thermal_Conductivity.f90"
#   include "Control_Mod/Physics/Scalars_Diffusivity.f90"
#   include "Control_Mod/Physics/Scalars_Diffusivities.f90"

    ! Multiphase flow
#   include "Control_Mod/Physics/Interface_Tracking.f90"
#   include "Control_Mod/Physics/Particle_Tracking.f90"
#   include "Control_Mod/Physics/Mass_Transfer_Model.f90"
#   include "Control_Mod/Physics/Lee_Model_Coefficients.f90"

    ! Heat transfer
#   include "Control_Mod/Physics/Heat_Transfer.f90"
#   include "Control_Mod/Physics/Buoyancy.f90"
#   include "Control_Mod/Physics/Reference_Temperature.f90"
#   include "Control_Mod/Physics/Saturation_Temperature.f90"
#   include "Control_Mod/Physics/Volume_Expansion_Coefficient.f90"

    ! Turbulence
#   include "Control_Mod/Physics/Hybrid_Les_Rans_Switch.f90"
#   include "Control_Mod/Physics/Monin_Obukov.f90"
#   include "Control_Mod/Physics/Roughness_Coefficient.f90"
#   include "Control_Mod/Physics/Rough_Walls.f90"
#   include "Control_Mod/Physics/Smagorinsky_Constant.f90"
#   include "Control_Mod/Physics/Wale_Constant.f90"
#   include "Control_Mod/Physics/Turbulence_Model.f90"
#   include "Control_Mod/Physics/Turbulent_Prandtl_Number.f90"
#   include "Control_Mod/Physics/Turbulent_Schmidt_Number.f90"
#   include "Control_Mod/Physics/Turbulent_Heat_Flux_Model.f90"
#   include "Control_Mod/Physics/Turbulent_Scalar_Flux_Model.f90"

    ! Other environmental conditions
#   include "Control_Mod/Physics/Bulk_Velocities.f90"
#   include "Control_Mod/Physics/Gravitational_Vector.f90"
#   include "Control_Mod/Physics/Pressure_Drops.f90"
#   include "Control_Mod/Physics/Point_For_Monitoring_Planes.f90"
#   include "Control_Mod/Physics/Potential_Initialization.f90"
#   include "Control_Mod/Physics/Volume_Flow_Rates.f90"

    ! Multiphase
#   include "Control_Mod/Physics/Max_Particles.f90"
#   include "Control_Mod/Physics/Number_Of_Phases.f90"
#   include "Control_Mod/Physics/Number_Of_Swarm_Sub_Steps.f90"
#   include "Control_Mod/Physics/Swarm_Coefficient_Of_Restitution.f90"
#   include "Control_Mod/Physics/Swarm_Density.f90"
#   include "Control_Mod/Physics/Swarm_Diameter.f90"
#   include "Control_Mod/Physics/Swarm_Subgrid_Scale_Model.f90"
#   include "Control_Mod/Physics/Phase_Capacities.f90"
#   include "Control_Mod/Physics/Phase_Conductivities.f90"
#   include "Control_Mod/Physics/Phase_Densities.f90"
#   include "Control_Mod/Physics/Phase_Viscosities.f90"
#   include "Control_Mod/Physics/Reference_Density.f90"
#   include "Control_Mod/Physics/Surface_Tension.f90"
#   include "Control_Mod/Physics/Track_Front.f90"
#   include "Control_Mod/Physics/Track_Surface.f90"

    ! Scalars (like species, for example)
#   include "Control_Mod/Physics/Number_Of_Scalars.f90"

    ! Statistics
#   include "Control_Mod/Physics/Starting_Time_Step_For_Turb_Statistics.f90"
#   include "Control_Mod/Physics/Starting_Time_Step_For_Swarm_Computation.f90"
#   include "Control_Mod/Physics/Starting_Time_Step_For_Swarm_Statistics.f90"

    ! Extrapolation to walls
#   include "Control_Mod/Physics/Extrapolate_Temperature_Exp.f90"

  end module

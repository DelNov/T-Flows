!==============================================================================!
  module Vof_Mod
!------------------------------------------------------------------------------!
!   Definition of variables used for all VOF multiphase modelling              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Front_Mod
  use Surf_Mod
  use Turb_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !               !
  !   Vof class   !
  !               !
  !---------------!
  type Vof_Type

    type(Grid_Type),  pointer :: pnt_grid  ! grid for which it is defined
    type(Field_Type), pointer :: pnt_flow  ! flow field for which it is defined
    type(Front_Type)          :: Front     ! pointer to Front (simple surface)
    type(Surf_Type)           :: surf      ! pointer to surface

    ! Volume fraction (colour function) and its smooth variant
    type(Var_Type) :: fun
    type(Var_Type) :: smooth

    ! Additional variables for calculation of vof function (fun)
    real, allocatable :: beta_f(:)
    real, allocatable :: beta_c(:)
    real, allocatable :: c_d(:)

    ! Surface curvature
    real, allocatable :: curv(:)

    ! Surface normals
    real, allocatable :: nx(:), ny(:), nz(:)

    ! Surface tension force
    real, allocatable :: surf_fx(:), surf_fy(:), surf_fz(:)  ! [N/m^3]

    ! Physical properties in case of multiphase flow
    real, allocatable :: phase_visc(:), phase_dens(:)
    real, allocatable :: phase_capa(:), phase_cond(:)
    real              :: surface_tension

    ! Phase change (called mass transfer to be consistent
    ! with heat transfer in the rest of the code)
    real    :: m_d, m_ini, m_s, m_s_acc

    ! Skewness correction
    logical :: skew_corr

    ! For phase change
    real :: t_sat, latent_heat  ! [K, J/kg]

    ! Heat from phase change and index of saturated cells
    real,    allocatable :: qci(:)
    real,    allocatable :: m_dot(:)         ! [kg/m^2 s]
    real                 :: add_mass_in, add_mass_out, vol_flux_avg
    type(Var_Type)       :: var

    ! User define parameters for vof function (fun)
    real    :: courant_max_param
    integer :: n_sub_param, corr_num_max
    integer :: n_conv_curv, n_conv_norm

    ! Averaging
    integer, allocatable :: avg_cells(:,:)

    ! Triangulate the front
    logical :: track_front

    ! Variable holding the multiphase model
    integer :: model

    contains

      !----------------------------------------!
      !   Procedures to advance vof function   !
      !----------------------------------------!
      procedure          :: Allocate_Vof
      procedure          :: Main_Vof
      procedure, private :: Compute_Vof
      procedure, private :: Discretize
      procedure, private :: Correct_Beta
      procedure, private :: Find_Upstream_Phi
      procedure, private :: Max_Courant_Number
      procedure, private :: Predict_Beta
      procedure, private :: Solve_System

      !------------------------------------------!
      !   Procedures helping to prepare smooth   !
      !   (convoluted) variant of vof function   !
      !   for eventual estimation of curvature   !
      !------------------------------------------!
      procedure, private :: Curvature_Csf
      procedure, private :: Smooth_Curvature
      procedure, private :: Smooth_For_Curvature_Csf
      procedure, private :: Smooth_Scalar

      !----------------------------------------------!
      !   Procedures to be called by other modules   !
      !----------------------------------------------!
      procedure :: Calculate_Grad_Matrix_With_Front
      procedure :: Get_Gas_And_Liquid_Phase
      procedure :: Grad_Component_No_Refresh_With_Front
      procedure :: Grad_Variable_With_Front
      procedure :: Mass_Transfer_Added_Volume
      procedure :: Mass_Transfer_Pressure_Source
      procedure :: Mass_Transfer_Vof_Source
      procedure :: Surface_Tension_Force
      procedure :: Update_Physical_Properties

  end type

  !--------------------------------------------------------!
  !   Parameters and variables defining multiphase model   !
  !--------------------------------------------------------!

  ! Parameters describing multiphase model choice
  ! (Prime numbers starting from 40000)
  integer, parameter :: NO_MULTIPHASE_MODEL   = 50021
  integer, parameter :: VOLUME_OF_FLUID       = 50023
  integer, parameter :: LAGRANGIAN_PARTICLES  = 50033
  integer, parameter :: EULER_EULER           = 50047
  integer, parameter :: FRONT_TRACKING        = 50051

  contains

  !----------------------------------------!
  !   Procedures to advance vof function   !
  !----------------------------------------!
  include 'Vof_Mod/Core/Allocate_Vof.f90'
  include 'Vof_Mod/Core/Main_Vof.f90'
  include 'Vof_Mod/Core/Compute_Vof.f90'
  include 'Vof_Mod/Core/Discretize.f90'
  include 'Vof_Mod/Core/Correct_Beta.f90'
  include 'Vof_Mod/Core/Find_Upstream_Phi.f90'
  include 'Vof_Mod/Core/Mass_Transfer_Added_Volume.f90'
  include 'Vof_Mod/Core/Max_Courant_Number.f90'
  include 'Vof_Mod/Core/Predict_Beta.f90'
  include 'Vof_Mod/Core/Solve_System.f90'

  !------------------------------------------!
  !   Procedures helping to prepare smooth   !
  !   (convoluted) variant of vof function   !
  !   for eventual estimation of curvature   !
  !------------------------------------------!
  include 'Vof_Mod/Curvature/Curvature_Csf.f90'
  include 'Vof_Mod/Curvature/Smooth_Curvature.f90'
  include 'Vof_Mod/Curvature/Smooth_For_Curvature_Csf.f90'
  include 'Vof_Mod/Curvature/Smooth_Scalar.f90'

  !----------------------------------------------!
  !   Procedures to be called by other modules   !
  !----------------------------------------------!
  include 'Vof_Mod/Utilities/Calculate_Grad_Matrix_With_Front.f90'
  include 'Vof_Mod/Utilities/Get_Gas_And_Liquid_Phase.f90'
  include 'Vof_Mod/Utilities/Grad_Component_No_Refresh_With_Front.f90'
  include 'Vof_Mod/Utilities/Grad_Variable_With_Front.f90'
  include 'Vof_Mod/Utilities/Mass_Transfer_Pressure_Source.f90'
  include 'Vof_Mod/Utilities/Mass_Transfer_Vof_Source.f90'
  include 'Vof_Mod/Utilities/Pressure_Correction.f90'
  include 'Vof_Mod/Utilities/Surface_Tension_Force.f90'
  include 'Vof_Mod/Utilities/Update_Physical_Properties.f90'

  end module

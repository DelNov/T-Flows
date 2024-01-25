#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Vof_Mod
!------------------------------------------------------------------------------!
!   Definition of variables used for all VOF multiphase modelling              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Stl_Mod
  use Surf_Mod  ! inherited from Front_Mod
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

    ! Stores the name of the STL file for initialization
    character(SL) :: name_stl = ''
    logical       :: init_stl = .false.     ! is it intialized from STL?

    type(Grid_Type),   pointer :: pnt_grid    ! its grid
    type(Field_Type),  pointer :: pnt_flow    ! its flow
    type(Matrix_Type), pointer :: pnt_matrix  ! its matrix
    type(Front_Type)           :: Front       ! pointer to front
    type(Surf_Type)            :: Surf        ! pointer to surface

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
    real :: phase_visc(0:1), phase_dens(0:1)
    real :: phase_capa(0:1), phase_cond(0:1)
    real :: surface_tension

    ! Phase change (called mass transfer to be consistent
    ! with heat transfer in the rest of the code)
    real    :: m_d, m_ini, m_s, m_s_acc

    ! Skewness correction
    logical :: skew_corr

    ! For phase change
    real :: t_sat, latent_heat  ! [K, J/kg]

    ! Neighbouring coefficients due to presence of saturated front
    real, allocatable :: a12(:)
    real, allocatable :: a21(:)
    real, allocatable :: m_dot(:)  ! [kg/s]

    type(Var_Type) :: t_0, t_1

    ! User define parameters for vof function (fun)
    real    :: courant_max_param
    integer :: n_sub_param, corr_num_max
    integer :: n_smooth_for_curv_csf

    ! Averaging
    integer, allocatable :: avg_cells(:,:)

    ! Mesh the front or surface
    logical :: track_front    ! simple but ugly
    logical :: track_surface  ! complicated but beautiful

    ! Coefficient of Lee for mass transfer rate
    ! m_dot [kg/s] = c_lee(1) * (T - Tsat) * A * V**(1/3) / latent_heat  for T <  Tsat
    !              = c_lee(2) * (T - Tsat) * A * V**(1/3) / latent_heat  for T >= Tsat
    ! m_dot [kg/s]: positive for vaporization, negative for condensation
    ! c_lee [1/s] : (1) coefficient for condensation, (2) for vaporisation
    ! A [m2]      : area of interface in cell
    ! V [m3]      : cell volume
    ! phase(0) and phase(1) are supporsed to be vapor and liquid phases, respectively.
    ! c_lee(1): phase(0) -> phase(1)
    ! c_lee(2): phase(1) -> phase(0)
    real    :: c_lee(2)       ! [1/s] 

    contains

      !----------------------------------------!
      !   Procedures to advance vof function   !
      !----------------------------------------!
      procedure          :: Create_Vof
      procedure          :: Main_Vof
      procedure, private :: Compute_Vof
      procedure, private :: Discretize
      procedure, private :: Correct_Beta
      procedure          :: Initialize_From_Stl
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
      procedure, private :: Smooth_Scalar
      procedure          :: Smooth_Vof_And_Compute_Surface_Normals

      !----------------------------------------------!
      !   Procedures to be called by other modules   !
      !----------------------------------------------!
      procedure :: Calculate_Grad_Matrix_With_Front
      procedure :: Extrapolate_Normal_To_Front
      procedure :: Get_Vapour_And_Liquid_Phase
      procedure :: Grad_Component_No_Refresh_With_Front
      procedure :: Grad_Variable_With_Front
      procedure :: Mass_Transfer_Added_Volume
      procedure :: Mass_Transfer_Estimate
      procedure :: Mass_Transfer_Pressure_Source
      procedure :: Mass_Transfer_Vof_Source
      procedure :: Surface_Tension_Force
      procedure :: Update_Physical_Properties

  end type

  contains

    !----------------------------------------!
    !   Procedures to advance vof function   !
    !----------------------------------------!
#   include "Vof_Mod/Core/Create_Vof.f90"
#   include "Vof_Mod/Core/Main_Vof.f90"
#   include "Vof_Mod/Core/Compute_Vof.f90"
#   include "Vof_Mod/Core/Discretize.f90"
#   include "Vof_Mod/Core/Correct_Beta.f90"
#   include "Vof_Mod/Core/Initialize_From_Stl.f90"
#   include "Vof_Mod/Core/Max_Courant_Number.f90"
#   include "Vof_Mod/Core/Predict_Beta.f90"
#   include "Vof_Mod/Core/Solve_System.f90"

    !------------------------------------------!
    !   Procedures helping to prepare smooth   !
    !   (convoluted) variant of vof function   !
    !   for eventual estimation of curvature   !
    !------------------------------------------!
#   include "Vof_Mod/Curvature/Curvature_Csf.f90"
#   include "Vof_Mod/Curvature/Smooth_Curvature.f90"
#   include "Vof_Mod/Curvature/Smooth_Vof_And_Compute_Surface_Normals.f90"
#   include "Vof_Mod/Curvature/Smooth_Scalar.f90"

    !----------------------------------------------!
    !   Procedures to be called by other modules   !
    !----------------------------------------------!
#   include "Vof_Mod/Utilities/Calculate_Grad_Matrix_With_Front.f90"
#   include "Vof_Mod/Utilities/Extrapolate_Normal_To_Front.f90"
#   include "Vof_Mod/Utilities/Get_Gas_And_Liquid_Phase.f90"
#   include "Vof_Mod/Utilities/Grad_Component_No_Refresh_With_Front.f90"
#   include "Vof_Mod/Utilities/Grad_Variable_With_Front.f90"
#   include "Vof_Mod/Utilities/Mass_Transfer_Added_Volume.f90"
#   include "Vof_Mod/Utilities/Mass_Transfer_Estimate.f90"
#   include "Vof_Mod/Utilities/Mass_Transfer_Pressure_Source.f90"
#   include "Vof_Mod/Utilities/Mass_Transfer_Vof_Source.f90"
#   include "Vof_Mod/Utilities/Surface_Tension_Force.f90"
#   include "Vof_Mod/Utilities/Update_Physical_Properties.f90"

    !------------------------------------!
    !   User functions for this module   !
    !------------------------------------!
#   include "User_Mod/Beginning_Of_Compute_Vof.f90"
#   include "User_Mod/End_Of_Compute_Vof.f90"

  end module

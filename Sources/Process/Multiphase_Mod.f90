!==============================================================================!
  module Multiphase_Mod
!------------------------------------------------------------------------------!
!   Definition of variables used for all multiphase modelling paradigms.       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Surf_Mod
  use Turb_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------------------!
  !                           !
  !      Multiphase type      !
  !   (variables needed for   !
  !   multiphase modelling)   !
  !                           !
  !---------------------------!
  type Multiphase_Type

    type(Grid_Type),  pointer :: pnt_grid  ! grid for which it is defined
    type(Field_Type), pointer :: pnt_flow  ! flow field for which it is defined
    type(Surf_Type)           :: surf      ! pointer to surface

    ! Volume fraction (colour function) and its smooth variant
    type(Var_Type) :: vof
    type(Var_Type) :: smooth

    ! Surface curvature
    real, allocatable :: curv(:)

    ! Surface normals
    real, allocatable :: nx(:), ny(:), nz(:)

    ! Surface tension force
    real, allocatable :: surf_fx(:), surf_fy(:), surf_fz(:)

    ! Physical properties in case of multiphase flow
    real, allocatable :: phase_visc(:), phase_dens(:)
    real, allocatable :: phase_capa(:), phase_cond(:)
    real              :: surface_tension

    ! Phase change
    real    :: m_d, m_ini, m_s, m_s_acc
    logical :: phase_change

    ! Skewness correction
    logical :: skew_corr

    ! For phase change
    real :: t_sat, latent_heat  ![K, J/kg]

    ! Heat from phase change and index of saturated cells
    real, allocatable    :: qci(:)
    real, allocatable    :: flux_rate(:)
    integer, allocatable :: ic(:)
    real                 :: add_mass_in, add_mass_out, vol_flux_avg

    ! User define parameters for vof
    real    :: courant_max_param
    integer :: n_sub_param, corr_num_max
    integer :: n_conv_curv, n_conv_norm

    ! Averaging
    integer, allocatable :: avg_cells(:,:)

    ! Triangulate the front
    logical :: track_front

    ! Variable holding the multiphase model
    integer :: model

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

  include 'Multiphase_Mod/Vof_Main.f90'
  include 'Multiphase_Mod/Vof_Allocate.f90'
  include 'Multiphase_Mod/Vof_Averaging.f90'
  include 'Multiphase_Mod/Vof_Compute.f90'
  include 'Multiphase_Mod/Vof_Coefficients.f90'
  include 'Multiphase_Mod/Vof_Correct_Beta.f90'
  include 'Multiphase_Mod/Vof_Curvature_Csf.f90'
  include 'Multiphase_Mod/Vof_Find_Upstream_Phi.f90'
  include 'Multiphase_Mod/Vof_Mass_Transfer.f90'
  include 'Multiphase_Mod/Vof_Mass_Transfer_Rate_In.f90'
  include 'Multiphase_Mod/Vof_Max_Courant_Number.f90'
  include 'Multiphase_Mod/Vof_Momentum_Contribution.f90'
  include 'Multiphase_Mod/Vof_Physical_Properties.f90'
  include 'Multiphase_Mod/Vof_Predict_Beta.f90'
  include 'Multiphase_Mod/Vof_Pressure_Correction.f90'
  include 'Multiphase_Mod/Vof_Smooth_Scalar.f90'
  include 'Multiphase_Mod/Vof_Smooth_Curvature.f90'
  include 'Multiphase_Mod/Vof_Solve_System.f90'
  include 'Multiphase_Mod/Vof_Surface_Tension_Contribution_Csf.f90'

  end module

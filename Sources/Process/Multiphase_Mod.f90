!==============================================================================!
  module Multiphase_Mod
!------------------------------------------------------------------------------!
!   Definition of variables used for all multiphase modelling paradigms.       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Var_Mod
  use Math_Mod
  use Face_Mod
  use Grid_Mod
  use Field_Mod
  use Cpu_Timer_Mod
  use Info_Mod
  use Solver_Mod
  use Control_Mod
  use Numerics_Mod
  use Const_Mod
  use Comm_Mod
  use Bulk_Mod
  use Matrix_Mod
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

    ! Volume fraction (colour function)
    type(Var_Type)    :: vof
    real, allocatable :: vof_f(:)
    real, allocatable :: curv(:)   ! curvature

    ! Distance function
    type(Var_Type)    :: dist_func

    ! Physical properties in case of multiphase flow
    real, allocatable :: phase_visc(:), phase_dens(:)
    real, allocatable :: phase_capa(:), phase_cond(:)
    real              :: surface_tension

    ! Phase change
    real              :: m_d, m_ini, m_s, m_s_acc
    logical           :: phase_change

    ! Skewness correction
    logical           :: skew_corr

    ! For calculation of distance function
    logical           :: d_func

    ! For phase change
    real              :: t_sat, latent_heat  ![K, J/kg]

    ! Surface force (at faces)
    real, allocatable :: fs_x(:)
    real, allocatable :: fs_y(:)
    real, allocatable :: fs_z(:)

    ! Surface force (at cells)
    real, allocatable :: fc_x(:)
    real, allocatable :: fc_y(:)
    real, allocatable :: fc_z(:)

    ! heat from phase change and index of saturated cells
    real, allocatable    :: qci(:)
    real, allocatable    :: flux_rate(:)
    integer, allocatable :: ic(:)
    real                 :: add_mass_in, add_mass_out, vol_flux_avg

    ! User define parameters for vof
    real    :: courant_max_param
    integer :: n_sub_param, corr_num_max
    integer :: n_conv_curv, n_conv_norm

    ! User defined parameters for distance function
    integer :: t_dist_scheme
    real    :: c_tau, c_eps

    ! Averaging
    integer, allocatable  :: avg_cells(:,:)

    ! Switch calculation curvature at nodes or at cells
    logical :: nodal_curvature

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

  contains

  include 'Multiphase_Mod/Alias_Vof.f90'
  include 'Multiphase_Mod/Allocate.f90'
  include 'Multiphase_Mod/Compute_Distance_Function.f90'
  include 'Multiphase_Mod/Compute_Vof.f90'
  include 'Multiphase_Mod/Vof_Averaging.f90'
  include 'Multiphase_Mod/Vof_Boundary_Extrapolation.f90'
  include 'Multiphase_Mod/Vof_Coefficients.f90'
  include 'Multiphase_Mod/Vof_Correct_Beta.f90'
  include 'Multiphase_Mod/Vof_Curvature_Csf.f90'
  include 'Multiphase_Mod/Vof_Curvature_Nodal.f90'
  include 'Multiphase_Mod/Vof_Find_Upstream_Phi.f90'
  include 'Multiphase_Mod/Vof_Find_Weight_Grad_From_Nodes.f90'
  include 'Multiphase_Mod/Vof_Find_Weight_Nodal_Grad.f90'
  include 'Multiphase_Mod/Vof_Grad_Component.f90'
  include 'Multiphase_Mod/Vof_Gradient_At_Nodes.f90'
  include 'Multiphase_Mod/Vof_Heavyside_Function.f90'
  include 'Multiphase_Mod/Vof_Mass_Transfer.f90'
  include 'Multiphase_Mod/Vof_Mass_Transfer_Rate_In.f90'
  include 'Multiphase_Mod/Vof_Max_Courant_Number.f90'
  include 'Multiphase_Mod/Vof_Momentum_Contribution.f90'
  include 'Multiphase_Mod/Vof_Nodal_Gradient.f90'
  include 'Multiphase_Mod/Vof_Open_Boundary.f90'
  include 'Multiphase_Mod/Vof_Predict_Beta.f90'
  include 'Multiphase_Mod/Vof_Pressure_Correction.f90'
  include 'Multiphase_Mod/Vof_Scale_Residuals.f90'
  include 'Multiphase_Mod/Vof_Skewness_Correction.f90'
  include 'Multiphase_Mod/Vof_Smooth_Scalar.f90'
  include 'Multiphase_Mod/Vof_Smooth_Curvature.f90'
  include 'Multiphase_Mod/Vof_Solve_System.f90'
  include 'Multiphase_Mod/Vof_Solver_Dist_Function_Cell_Loop.f90'
  include 'Multiphase_Mod/Vof_Surface_Tension_Contribution.f90'
  include 'Multiphase_Mod/Vof_Surface_Tension_Contribution_Nodal.f90'
  include 'Multiphase_Mod/Update_Physical_Properties.f90'

  end module

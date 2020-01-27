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
    real              :: surface_tension

    ! For calculation of distance function
    logical           :: d_func

    ! Body force
    real, allocatable :: body_fx(:)
    real, allocatable :: body_fy(:)
    real, allocatable :: body_fz(:)

    ! User define parameters for vof
    real              :: courant_max_param
    integer           :: n_sub_param, corr_num_max
    integer           :: n_conv_curv, n_conv_norm

    ! User defined parameters for distance function
    integer           :: t_dist_scheme
    real              :: c_tau, c_eps
  end type

  !--------------------------------------------------------!
  !   Parameters and variables defining multiphase model   !
  !--------------------------------------------------------!

  ! Variable holding the multiphase model
  integer :: multiphase_model

  ! Parameters describing multiphase model choice
  ! (Prime numbers starting from 40000)
  integer, parameter :: NONE                  = 50021
  integer, parameter :: VOLUME_OF_FLUID       = 50023
  integer, parameter :: LAGRANGIAN_PARTICLES  = 50033
  integer, parameter :: EULER_EULER           = 50047

  contains

  include 'Multiphase_Mod/Alias_Vof.f90'
  include 'Multiphase_Mod/Allocate.f90'
  include 'Multiphase_Mod/Compute_Distance_Function.f90'
  include 'Multiphase_Mod/Compute_Vof.f90'
  include 'Multiphase_Mod/Vof_Coefficients.f90'
  include 'Multiphase_Mod/Vof_Correct_Beta.f90'
  include 'Multiphase_Mod/Vof_Find_Upstream_Phi.f90'
  include 'Multiphase_Mod/Vof_Heavyside_Function.f90'
  include 'Multiphase_Mod/Vof_Predict_Beta.f90'
  include 'Multiphase_Mod/Vof_Pressure_Correction.f90'
  include 'Multiphase_Mod/Vof_Max_Courant_Number.f90'
  include 'Multiphase_Mod/Vof_Momentum_Contribution.f90'
  include 'Multiphase_Mod/Vof_Solve_System.f90'
  include 'Multiphase_Mod/Vof_Solver_Dist_Function_Cell_Loop.f90'
  include 'Multiphase_Mod/Vof_Surface_Tension_Contribution.f90'
  include 'Multiphase_Mod/Update_Physical_Properties.f90'

  end module

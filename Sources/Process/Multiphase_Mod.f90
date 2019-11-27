!==============================================================================!
  module Multiphase_Mod
!------------------------------------------------------------------------------!
!   Definition of variables used for all multiphase modelling paradigms.       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Var_Mod,       only: Var_Type, Var_Mod_Allocate_Solution
  use Math_Mod
  use Face_Mod,      only: Face_Type
  use Grid_Mod,      only: Grid_Type
  use Field_Mod,     only: Field_Type, density, viscosity, dens_face
  use Cpu_Timer_Mod, only: Cpu_Timer_Mod_Start, Cpu_Timer_Mod_Stop
  use Info_Mod,      only: Info_Mod_Iter_Fill_At
  use Solver_Mod,    only: Solver_Type, Bicg, Cg, Cgs, Acm
  use Grad_Mod
  use Control_Mod
  use Numerics_Mod
  use Const_Mod
  use Comm_Mod
  use Bulk_Mod,      only: Bulk_Type
  use Matrix_Mod,    only: Matrix_Type
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

  end type

  ! Physical properties in case of multiphase flow
  real, allocatable :: phase_visc(:), phase_dens(:)
  real              :: surface_tension

  ! Body force
  real, allocatable :: body_fx(:)
  real, allocatable :: body_fy(:)
  real, allocatable :: body_fz(:)


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
  include 'Multiphase_Mod/Compute_Vof.f90'
  include 'Multiphase_Mod/Update_Physical_Properties.f90'
  include 'Multiphase_Mod/Vof_Correct_Beta.f90'
  include 'Multiphase_Mod/Vof_Predict_Beta.f90'
  include 'Multiphase_Mod/Vof_Spurious_Post.f90'
  include 'Multiphase_Mod/Vof_Surface_Tension_Contribution.f90'
  include 'Multiphase_Mod/Vof_Max_Courant_Number.f90'
  include 'Multiphase_Mod/Vof_Pressure_Correction.f90'
  include 'Multiphase_Mod/Vof_Momentum_Contribution.f90'
  include 'Multiphase_Mod/Vof_Coefficients.f90'
  include 'Multiphase_Mod/Vof_Solve_System.f90'

  end module

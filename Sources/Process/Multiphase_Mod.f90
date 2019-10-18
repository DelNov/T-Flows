!==============================================================================!
  module Multiphase_Mod
!------------------------------------------------------------------------------!
!   Definition of variables used for all multiphase modelling paradigms.       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Var_Mod,   only: Var_Type, Var_Mod_Allocate_Solution
  use Grid_Mod,  only: Grid_Type
  use Field_Mod, only: Field_Type, density, viscosity, dens_face
  use Grad_Mod
  use Control_Mod
  use Numerics_Mod
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

  end type

  ! Physical properties in case of multiphase flow
  real, allocatable :: phase_visc(:), phase_dens(:)
  real              :: surface_tension

  !--------------------------------------------------------!
  !   Parameters and variables defining multiphase model   !
  !--------------------------------------------------------!

  ! Variable holding the multiphase model
  integer :: multiphase_model

  ! Parameters describing multiphase model choice
  ! (Prime numbers starting from 40000)
  integer, parameter :: VOLUME_OF_FLUID       = 40009
  integer, parameter :: LAGRANGIAN_PARTICLES  = 40013
  integer, parameter :: EULER_EULER           = 40031

  contains

  include 'Multiphase_Mod/Alias_Vof.f90'
  include 'Multiphase_Mod/Allocate.f90'
  include 'Multiphase_Mod/Compute_Vof.f90'
  include 'Multiphase_Mod/Vof_Correct_Beta.f90'
  include 'Multiphase_Mod/Vof_Initialization.f90'
  include 'Multiphase_Mod/Vof_Initialization_Cylinder.f90'
  include 'Multiphase_Mod/Vof_Initialization_Ellipsoid.f90'
  include 'Multiphase_Mod/Vof_Initialization_Plane.f90'
  include 'Multiphase_Mod/Vof_Predict_Beta.f90'
  include 'Multiphase_Mod/Vof_Spurious_Post.f90'
  include 'Multiphase_Mod/Vof_Surface_Tension_Contribution.f90'

  end module

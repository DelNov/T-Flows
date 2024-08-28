#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"
#include "../Shared/Macros.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Turb_Mod
!------------------------------------------------------------------------------!
!>  Definition of variables used for all turbulence modelling paradigms.       !
!>
!>  This is s a simplified version from the same module in Process_Cpu.
!>  Hopefully, as more modules are ported to Process_Gpu, this module will
!>  get closer and closer to its origin from Process_Cpu.  In any case, it
!>  shoiuld be our utmost priority to keep this file and its counterpart
!>  in Process_Cpu as close to one another as possible.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Info_Mod
  use Iter_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------------------!
  !                           !
  !      Turbulence type      !
  !   (variables needed for   !
  !   turbulence modelling)   !
  !                           !
  !---------------------------!
  type Turb_Type

    !---------------!
    !   Variables   !
    !---------------!

    ! Non-dimensional distance
    real, allocatable :: y_plus(:)

    ! Turbulent viscosity
    real, allocatable :: vis_t(:) ! [kg/(m s)]

    ! Scale-resolving simulations
    real, allocatable :: wale_v(:)

    ! Variable holding the turbulence model; its variant and statistics
    integer :: model

    contains
      procedure :: Init_Turb
      procedure :: Main_Turb
      procedure :: Create_Turb

      ! Computation of turbulence viscosity
      procedure, private :: Vis_T_Subgrid
      procedure, private :: Vis_T_Wale

  end type

  ! Parameters describing turbulence model choice
  ! (Prime numbers starting from 30000)
  integer, parameter :: NO_TURBULENCE_MODEL   = 30011
  integer, parameter :: LES_SMAGORINSKY       = 30029
  integer, parameter :: LES_WALE              = 30059

  !--------------------------------!
  !   Turbulence model constants   !
  !--------------------------------!

  ! For scale-resolving models
  real :: c_smag

  !-----------------------------------!
  !   Auxiliary turbulent variables   !
  !-----------------------------------!

  contains

    ! Logic of turbulence models
#   include "Turb_Mod/Init_Turb.f90"
#   include "Turb_Mod/Main_Turb.f90"

    ! The constructor-like
#   include "Turb_Mod/Create_Turb.f90"


    ! Functions to set turbulence constants
#   include "Turb_Mod/Const_Les.f90"
    ! Computation of turbulence viscosity
#   include "Turb_Mod/Vis_T_Subgrid.f90"
#   include "Turb_Mod/Vis_T_Wale.f90"

  end module

!==============================================================================!
  module Numerics_Mod
!------------------------------------------------------------------------------!
!   Module for basic flow field plus temperature.                              !
!   It is a bit of a mumbo-jumbo at this moment, it will furhter have to       !
!   differentiate into numerical and physica parts.                            !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! These paramters should not be here but in a new module Numerics_Mod
  integer, parameter :: UPWIND    = 40009
  integer, parameter :: CENTRAL   = 40013
  integer, parameter :: LUDS      = 40031
  integer, parameter :: QUICK     = 40037
  integer, parameter :: SMART     = 40039
  integer, parameter :: GAMMA     = 40063
  integer, parameter :: MINMOD    = 40087
  integer, parameter :: BLENDED   = 40093
  integer, parameter :: SUPERBEE  = 40099   
  integer, parameter :: AVL_SMART = 40111

  ! Variable holding time integration scheme
  integer :: time_innertial
  integer :: time_advection
  integer :: time_diffusion
  integer :: time_cross_diff

  ! Time integration parameters
  integer, parameter :: FULLY_IMPLICIT  = 50021
  integer, parameter :: ADAMS_BASHFORTH = 50023
  integer, parameter :: CRANK_NICOLSON  = 50033
  integer, parameter :: LINEAR          = 50047
  integer, parameter :: PARABOLIC       = 50051

  end module

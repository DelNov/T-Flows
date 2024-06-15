#include "../Shared/Browse.h90"

!==============================================================================!
  module Numerics_Mod
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Message_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Time integration parameters
  integer, parameter :: LINEAR        = 40129
  integer, parameter :: PARABOLIC     = 40151

  contains

#   include "Numerics_Mod/Time_Integration_Scheme_Code.f90"

  end module

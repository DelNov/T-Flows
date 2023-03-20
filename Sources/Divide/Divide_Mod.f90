#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Divide_Mod
!------------------------------------------------------------------------------!
!   Collection of functions used in the Divide program.  In honesty, it was    !
!   introduced to get rid of the Fortran header files with interfaces which,   !
!   in effect, was needed for Intel compiler to work in the debug mode.        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------!
  !   Divide type   !
  !-----------------!
  type Divide_Type

    contains
      procedure :: Save_Subdomains
      procedure :: Logo_Div
  end type

  type(Divide_Type) :: Divide

  contains

#   include "Divide_Mod/Save_Subdomains.f90"
#   include "Divide_Mod/Logo_Div.f90"

  end module

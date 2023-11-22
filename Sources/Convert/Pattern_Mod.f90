#include "../Shared/Assert.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Pattern_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !   Pattern type   !
  !------------------!
  type Pattern_Type
    integer(1), allocatable :: pattern(:)
    integer                 :: length
    contains
      procedure :: Create_Pattern
      procedure :: Match_Pattern
  end type

  contains

#   include "Pattern_Mod/Create_Pattern.f90"
#   include "Pattern_Mod/Match_Pattern.f90"

  end module

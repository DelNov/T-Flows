!==============================================================================!
  module String_Mod
!------------------------------------------------------------------------------!
!   Procedures for handling string.  For now To_Upper_Case and To_Lower_Case   !
!   and a function which returns a string with only first letter in upper.     !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------!
  !   String type   !
  !-----------------!
  type String_Type

    contains
      procedure :: First_Upper
      procedure :: To_Lower_Case
      procedure :: To_Upper_Case

  end type

  ! One handle for all modules to use string routines
  type(String_Type) :: String

  contains

#   include "String_Mod/First_Upper.f90"
#   include "String_Mod/To_Lower_Case.f90"
#   include "String_Mod/To_Upper_Case.f90"

  end module

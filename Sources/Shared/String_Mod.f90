#include "../Shared/Unused.h90"

!==============================================================================!
  module String_Mod
!------------------------------------------------------------------------------!
!>  The String_Mod module is a simple and straightforward module designed for
!>  handling string transformations. It includes procedures for converting
!>  strings to upper case, lower case, and a function to capitalize only the
!>  first letter of a string.
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------!
  !   String type   !
  !-----------------!
  !> String_Type encapsulates procedures related to string manipulation.
  type String_Type

    contains
      procedure :: First_Upper
      procedure :: To_Lower_Case
      procedure :: To_Upper_Case

  end type

  ! One handle for all modules to use string routines
  type(String_Type) :: String  !! A global object of String_Type that can be
                               !! used across different parts of the program.
  contains

#   include "String_Mod/First_Upper.f90"
#   include "String_Mod/To_Lower_Case.f90"
#   include "String_Mod/To_Upper_Case.f90"

  end module

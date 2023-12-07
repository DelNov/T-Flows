#include "../Shared/Assert.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Memory_Mod
!------------------------------------------------------------------------------!
!>  This module provides a robust framework for dynamic memory management in
!>  Fortran. It encapsulates methods for safely handling, resizing, and probing
!>  both arrays and matrices, supporting both integer and real data types. It
!>  is designed around a Memory_Type type that contains all the necessary
!>  procedures for memory operations.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Message_Mod
  use Assert_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------!
  !   Memory type   !
  !-----------------!
  !> A custom type encapsulating memory management procedures.
  type Memory_Type

    contains
      procedure          :: Array_Int
      procedure          :: Array_Log
      procedure          :: Array_Real
      procedure          :: Matrix_Int
      procedure          :: Matrix_Log
      procedure          :: Matrix_Real
      procedure, private :: Test_Array_Int
      procedure, private :: Test_Array_Log
      procedure, private :: Test_Array_Real
      procedure, private :: Test_Matrix_Int
      procedure, private :: Test_Matrix_Log
      procedure, private :: Test_Matrix_Real
      procedure, private :: Work_Out_I_Ranges
      procedure, private :: Work_Out_J_Ranges

  end type

  type(Memory_Type) :: Enlarge  !! a singleton object Memory_Type used for
                                !! accessing memory management procedures

  contains
#   include "Memory_Mod/Array_Int.f90"
#   include "Memory_Mod/Array_Log.f90"
#   include "Memory_Mod/Array_Real.f90"
#   include "Memory_Mod/Matrix_Int.f90"
#   include "Memory_Mod/Matrix_Log.f90"
#   include "Memory_Mod/Matrix_Real.f90"
#   include "Memory_Mod/Test_Array_Int.f90"
#   include "Memory_Mod/Test_Array_Log.f90"
#   include "Memory_Mod/Test_Array_Real.f90"
#   include "Memory_Mod/Test_Matrix_Int.f90"
#   include "Memory_Mod/Test_Matrix_Log.f90"
#   include "Memory_Mod/Test_Matrix_Real.f90"
#   include "Memory_Mod/Work_Out_I_Ranges.f90"
#   include "Memory_Mod/Work_Out_J_Ranges.f90"

  end module





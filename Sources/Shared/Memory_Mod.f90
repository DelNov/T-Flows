!==============================================================================!
  module Memory_Mod
!------------------------------------------------------------------------------!
!>  This module provides a robust framework for dynamic memory management in
!>  Fortran. It encapsulates methods for safely handling, resizing, and probing
!>  both arrays and matrices, supporting both integer and real data types. It
!>  is designed around a Memory_Type type that contains all the necessary
!>  procedures for memory operations.
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
      procedure          :: Array_Real
      procedure          :: Matrix_Int
      procedure          :: Matrix_Real
      procedure, private :: Probe_Int_Array
      procedure, private :: Probe_Int_Matrix
      procedure, private :: Probe_Real_Array
      procedure, private :: Probe_Real_Matrix

  end type

  type(Memory_Type) :: Enlarge  !! a singleton object Memory_Type used for
                                !! accessing memory management procedures

  contains
    include "Memory_Mod/Array_Int.f90"
    include "Memory_Mod/Array_Real.f90"
    include "Memory_Mod/Matrix_Int.f90"
    include "Memory_Mod/Matrix_Real.f90"
    include "Memory_Mod/Probe_Int_Array.f90"
    include "Memory_Mod/Probe_Int_Matrix.f90"
    include "Memory_Mod/Probe_Real_Array.f90"
    include "Memory_Mod/Probe_Real_Matrix.f90"

  end module





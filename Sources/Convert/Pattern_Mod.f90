#include "../Shared/Assert.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Pattern_Mod
!------------------------------------------------------------------------------!
!>  The Pattern_Mod module is designed for pattern recognition and matching in
!>  input files within the Convert sub-program. It is particularly useful in
!>  procedures Load_Gmsh, where identifying specific patterns in the input
!>  files is crucial for efficient reading and file processing.
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !   Pattern type   !
  !------------------!
  !> Encapsulates data and mamber procedures to facilitate
  !> search for patterns in input files.  Used in Load_Gmsh.
  type Pattern_Type
    integer(1), allocatable :: pattern(:)  !! storage for pattern
    integer                 :: length      !! pattern length
    contains
      procedure :: Create_Pattern
      procedure :: Match_Pattern
  end type

  contains

#   include "Pattern_Mod/Create_Pattern.f90"
#   include "Pattern_Mod/Match_Pattern.f90"

  end module

# include "../Shared/Unused.h90"

# define Approx(x, y)  Math % Approx_Real(x, y)

!==============================================================================!
  module Sort_Mod
!------------------------------------------------------------------------------!
!>  A holder of the Sort_Type, which holds a collection of sorting routines.
!>  If you need sort an array, have a look here first, a routine you need might
!>  already be defined here.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Math_Mod
  use Swap_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Sort type   !
  !---------------!
  !> Sort_Type holds a collection of routines to sort arrays of different types
  !> It comes in two flavors: heap-sort routines and quick-sort routines, which
  !> can be set through compiler option SORT=heap/quick.  Quick is a bit faster,
  !> but sometimes runs out of heap memory because of its recursive nature.  If
  !> that happens, recompile the programs with option SORT=heap.
  type Sort_Type

    contains
      procedure :: Int_Array
      procedure :: Int_By_Index
      procedure :: Int_Carry_Int
      procedure :: Int_Carry_Real
      procedure :: Real_Array
      procedure :: Real_By_Index
      procedure :: Real_Carry_Int
      procedure :: Real_Carry_Two_Int
      procedure :: Reverse_Order_Int
      procedure :: Reverse_Order_Real
      procedure :: Three_Int_Carry_Two_Int
      procedure :: Three_Int_Carry_Three_Int
      procedure :: Three_Int_Carry_Int
      procedure :: Three_Int                   ! unused
      procedure :: Three_Real_Carry_Two_Int
      procedure :: Three_Real_Carry_Three_Int
      procedure :: Three_Real_Carry_Int
      procedure :: Three_Real                  ! unused
      procedure :: Two_Int_Carry_Two_Int
      procedure :: Two_Int_Carry_Int
      procedure :: Two_Int                     ! unused
      procedure :: Two_Real_Carry_Two_Int
      procedure :: Two_Real_Carry_Int
      procedure :: Two_Real                    ! unused
      procedure :: Unique_Int

  end type

  !--------------------------------------!
  !   Create one instance of Sort type   !
  !     for all other modules to use     !
  !--------------------------------------!
  type(Sort_Type) :: Sort

  contains

#   include "Sort_Mod/Int_By_Index.f90"
#   include "Sort_Mod/Real_By_Index.f90"
#   include "Sort_Mod/Reverse_Order_Int.f90"
#   include "Sort_Mod/Reverse_Order_Real.f90"
#   include "Sort_Mod/Unique_Int.f90"

#   if T_FLOWS_QUICKSORT == 1
#     include "Sort_Mod/Quick/Int_Array.f90"
#     include "Sort_Mod/Quick/Int_Carry_Int.f90"
#     include "Sort_Mod/Quick/Int_Carry_Real.f90"
#     include "Sort_Mod/Quick/Real_Array.f90"
#     include "Sort_Mod/Quick/Real_Carry_Int.f90"
#     include "Sort_Mod/Quick/Real_Carry_Two_Int.f90"
#     include "Sort_Mod/Quick/Three_Int_Carry_Two_Int.f90"
#     include "Sort_Mod/Quick/Three_Int_Carry_Three_Int.f90"
#     include "Sort_Mod/Quick/Three_Int_Carry_Int.f90"
#     include "Sort_Mod/Quick/Three_Int.f90"
#     include "Sort_Mod/Quick/Three_Real_Carry_Two_Int.f90"
#     include "Sort_Mod/Quick/Three_Real_Carry_Three_Int.f90"
#     include "Sort_Mod/Quick/Three_Real_Carry_Int.f90"
#     include "Sort_Mod/Quick/Three_Real.f90"
#     include "Sort_Mod/Quick/Two_Int_Carry_Two_Int.f90"
#     include "Sort_Mod/Quick/Two_Int_Carry_Int.f90"
#     include "Sort_Mod/Quick/Two_Int.f90"
#     include "Sort_Mod/Quick/Two_Real_Carry_Two_Int.f90"
#     include "Sort_Mod/Quick/Two_Real_Carry_Int.f90"
#     include "Sort_Mod/Quick/Two_Real.f90"
#   else
#     include "Sort_Mod/Heap/Int_Array.f90"
#     include "Sort_Mod/Heap/Int_Carry_Int.f90"
#     include "Sort_Mod/Heap/Int_Carry_Real.f90"
#     include "Sort_Mod/Heap/Real_Array.f90"
#     include "Sort_Mod/Heap/Real_Carry_Int.f90"
#     include "Sort_Mod/Heap/Real_Carry_Two_Int.f90"
#     include "Sort_Mod/Heap/Three_Int_Carry_Two_Int.f90"
#     include "Sort_Mod/Heap/Three_Int_Carry_Three_Int.f90"
#     include "Sort_Mod/Heap/Three_Int_Carry_Int.f90"
#     include "Sort_Mod/Heap/Three_Int.f90"
#     include "Sort_Mod/Heap/Three_Real_Carry_Two_Int.f90"
#     include "Sort_Mod/Heap/Three_Real_Carry_Three_Int.f90"
#     include "Sort_Mod/Heap/Three_Real_Carry_Int.f90"
#     include "Sort_Mod/Heap/Three_Real.f90"
#     include "Sort_Mod/Heap/Two_Int_Carry_Two_Int.f90"
#     include "Sort_Mod/Heap/Two_Int_Carry_Int.f90"
#     include "Sort_Mod/Heap/Two_Int.f90"
#     include "Sort_Mod/Heap/Two_Real_Carry_Two_Int.f90"
#     include "Sort_Mod/Heap/Two_Real_Carry_Int.f90"
#     include "Sort_Mod/Heap/Two_Real.f90"
#   endif

  end module

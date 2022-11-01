!==============================================================================!
  module Sort_Mod
!------------------------------------------------------------------------------!
!   A collection of sorting (and maybe compression in the future) subroutines  !
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
      procedure :: Three_Int
      procedure :: Three_Real_Carry_Two_Int
      procedure :: Three_Real_Carry_Three_Int
      procedure :: Three_Real_Carry_Int
      procedure :: Three_Real
      procedure :: Two_Int_Carry_Two_Int
      procedure :: Two_Int_Carry_Int
      procedure :: Two_Int
      procedure :: Two_Real_Carry_Two_Int
      procedure :: Two_Real_Carry_Int
      procedure :: Two_Real
      procedure :: Unique_Int

  end type

  !--------------------------------------!
  !   Create one instance of Sort type   !
  !     for all other modules to use     !
  !--------------------------------------!
  type(Sort_Type) :: Sort

  contains

#   include "Sort_Mod/Int_Array.f90"
#   include "Sort_Mod/Int_By_Index.f90"
#   include "Sort_Mod/Int_Carry_Int.f90"
#   include "Sort_Mod/Int_Carry_Real.f90"
#   include "Sort_Mod/Real_Array.f90"
#   include "Sort_Mod/Real_By_Index.f90"
#   include "Sort_Mod/Real_Carry_Int.f90"
#   include "Sort_Mod/Real_Carry_Two_Int.f90"
#   include "Sort_Mod/Reverse_Order_Int.f90"
#   include "Sort_Mod/Reverse_Order_Real.f90"
#   include "Sort_Mod/Three_Int_Carry_Two_Int.f90"
#   include "Sort_Mod/Three_Int_Carry_Three_Int.f90"
#   include "Sort_Mod/Three_Int_Carry_Int.f90"
#   include "Sort_Mod/Three_Int.f90"
#   include "Sort_Mod/Three_Real_Carry_Two_Int.f90"
#   include "Sort_Mod/Three_Real_Carry_Three_Int.f90"
#   include "Sort_Mod/Three_Real_Carry_Int.f90"
#   include "Sort_Mod/Three_Real.f90"
#   include "Sort_Mod/Two_Int_Carry_Two_Int.f90"
#   include "Sort_Mod/Two_Int_Carry_Int.f90"
#   include "Sort_Mod/Two_Int.f90"
#   include "Sort_Mod/Two_Real_Carry_Two_Int.f90"
#   include "Sort_Mod/Two_Real_Carry_Int.f90"
#   include "Sort_Mod/Two_Real.f90"
#   include "Sort_Mod/Unique_Int.f90"

  end module

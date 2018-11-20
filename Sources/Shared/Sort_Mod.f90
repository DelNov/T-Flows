!==============================================================================!
  module Sort_Mod
!------------------------------------------------------------------------------!
!   A collection of sorting (and maybe compression) subroutines                !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  contains

  include 'Sort_Mod/2_Int_Carry_Int.f90'
  include 'Sort_Mod/2_Int.f90'
  include 'Sort_Mod/2_Real_Carry_Int.f90'
  include 'Sort_Mod/2_Real.f90'
  include 'Sort_Mod/3_Int_Carry_Int.f90'
  include 'Sort_Mod/3_Int.f90'
  include 'Sort_Mod/3_Real_Carry_Int.f90'
  include 'Sort_Mod/3_Real.f90'
  include 'Sort_Mod/Int_By_Index.f90'
  include 'Sort_Mod/Int_Carry_Int.f90'
  include 'Sort_Mod/Int_Carry_Real.f90'
  include 'Sort_Mod/Int.f90'
  include 'Sort_Mod/Real_By_Index.f90'
  include 'Sort_Mod/Real_Carry_Int.f90'
  include 'Sort_Mod/Real.f90'
  include 'Sort_Mod/Short_Carry_Short.f90'

  end module

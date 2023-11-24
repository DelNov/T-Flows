!==============================================================================!
  module Swap_Mod
!------------------------------------------------------------------------------!
!>  Swap_Mod is a utility module in T-Flows designed to hold subroutine for
!>  swapping values of two variables. It includes subroutines for swapping
!>  integers and real numbers.
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  contains

#   include "Swap_Mod/Swap_Int.f90"
#   include "Swap_Mod/Swap_Real.f90"

  end module

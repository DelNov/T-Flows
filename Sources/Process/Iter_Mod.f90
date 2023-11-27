#include "../Shared/Assert.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Iter_Mod
!------------------------------------------------------------------------------!
!>  This module encapsulates the data and procedures relevant to the iteration
!>  process in SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) and
!>  PISO (Pressure-Implicit with Splitting of Operators) algorithms used in
!>  T-Flows.  Module entials the definition of the Iter_Type, global singleton
!>  object Iter and a lot of iteration-specific data and methods.
!>  The Iter object is primarily used in the main function of the Process
!>  program to manage outer iterations within the SIMPLE or PISO algorithms.
!>  Due to its global accessibility, it can be utilized in other functions
!>  where decisions depend on the current iteration state.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Iter type   !
  !---------------!
  !> A type that holds iteration-specific data and methods and is
  !> the basis for definition of the global singleton object Iter
  type Iter_Type

    integer, private :: current_iteration = 0  !! current iteration
    integer, private :: max_iterations         !! maximum number of iterations
    integer, private :: min_iterations         !! minimum number of iterations
    real,    private :: tolerance              !! tolerance for iterations

    contains
      procedure :: Current
      procedure :: Get_Max
      procedure :: Get_Min
      procedure :: Get_Tol
      procedure :: Max_Fields_Residual
      procedure :: Needs_More_Iterations
      procedure :: Set_Max
      procedure :: Set_Min
      procedure :: Set_Tol

  end type

  !----------------------!
  !   Singleton object   !
  !----------------------!
  type(Time_Type) :: Iter  !! singleton type object of type Iter_Type defined
                           !! globally for synchronization of its data members
                           !! and easier access to its member functions
  contains

#   include "Iter_Mod/Current.f90"
#   include "Iter_Mod/Get_Max.f90"
#   include "Iter_Mod/Get_Min.f90"
#   include "Iter_Mod/Get_Tol.f90"
#   include "Iter_Mod/Max_Fields_Residual.f90"
#   include "Iter_Mod/Needs_More_Iterations.f90"
#   include "Iter_Mod/Set_Max.f90"
#   include "Iter_Mod/Set_Min.f90"
#   include "Iter_Mod/Set_Tol.f90"

  end module

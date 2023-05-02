#include "../Shared/Assert.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Iter_Mod
!------------------------------------------------------------------------------!
!   Module containing data and procedures pertinent to SIMPLE/PISO iterations  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Iter type   !
  !---------------!
  type Iter_Type

    integer, private :: current_iteration = 0
    integer, private :: max_iterations         ! maximum number of iterations
    integer, private :: min_iterations         ! minimum number of iterations
    real,    private :: tolerance              ! tolerance for iterations

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
  type(Iter_Type) :: Iter

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

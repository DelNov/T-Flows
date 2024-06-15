#include "../Shared/Assert.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Var_Mod
!------------------------------------------------------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

  !--------------!
  !   Var type   !
  !--------------!
  type Var_Type

    type(Grid_Type), pointer :: pnt_grid  !! the grid on which it's defined

    character(VL) :: name       !! variable name, always upper case and
                                !! very short (4, defined in Const_Mod)
    character(VL) :: flux_name  !! variable flux name, always upper case and
                                !! very short (4, defined in Const_Mod)
    real, allocatable :: n(:)   !! new value (current time step)
    real, allocatable :: o(:)   !! old value (previous time step)
    real, allocatable :: b(:)   !! boundary value
    real, allocatable :: q(:)   !! wall flux of the variable
    real              :: res    !! residual after linear solver

    ! Boundary cell type (important for scalars, since they
    ! can have different boundary conditions at the walls)
    integer, allocatable :: bnd_cond_type(:)

    ! Parameters for numerical solution of the variable
    real          :: blend         !! upwind blending (1.0 central; 0.0 upwind)
    real          :: tol = PICO    !! linear solver tolerance
    real          :: urf           !! under-relaxation factor
    integer       :: miter         !! max number of iterations for variable

  end type

  contains

#   include "Var_Mod/Create_Variable.f90"

  end module

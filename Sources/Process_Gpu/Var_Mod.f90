#include "../Shared/Assert.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Var_Mod
!------------------------------------------------------------------------------!
!>  Module Var_Mod defines Var_Type which is, in turn, used to define all
!>  dependent variables in the T-Flows' solver Process.  Thus, individual
!>  velocity components, temperature, pressure and pressure corrections are all
!>  objects of Var_Type.  The Var_Type holds data fields to store its current
!>  value (n for new) old value (o for old) and older than old (oo) value.
!>  In addition, it also holds many linear solver parameters, such as desired
!>  solver tolerance levels required in linear solvers.
!------------------------------------------------------------------------------!
  use Numerics_Mod
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
    real, allocatable :: oo(:)  !! time step before old (older than old)
    real, allocatable :: b(:)   !! boundary value
    real, allocatable :: q(:)   !! wall flux of the variable
    real              :: res    !! residual after linear solver

    ! Boundary cell type (important for scalars, since they
    ! can have different boundary conditions at the walls)
    integer, allocatable :: bnd_cond_type(:)

    ! Parameters for numerical solution of the variable
    real          :: blend         !! upwind blending (1.0 central; 0.0 upwind)
    integer       :: td_scheme     !! time-disretization scheme
    real          :: tol = PICO    !! linear solver tolerance
    real          :: urf           !! under-relaxation factor
    integer       :: miter         !! max number of iterations for variable
    integer       :: niter         !! executed number of iterations for var
    logical       :: blend_matrix  !! are you blending matrix with upwind?

  end type

  contains

#   include "Var_Mod/Create_Variable.f90"

  end module

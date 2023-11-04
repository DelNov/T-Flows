#include "../Shared/Assert.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Var_Mod
!------------------------------------------------------------------------------!
  use Petsc_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

  !--------------!
  !   Var type   !
  !--------------!
  type Var_Type

    type(Grid_Type),   pointer :: pnt_grid    ! its grid
    type(Matrix_Type), pointer :: pnt_matrix  ! its matrix

    character(VL) :: name       ! variable name, always upper case and
                                ! very short (4, defined in Const_Mod)
    character(VL) :: flux_name  ! variable flux name, always upper case and
                                ! very short (4, defined in Const_Mod)
    real, allocatable :: n(:)              ! new value
    real, allocatable :: o(:), oo(:)       ! old and older then old
    real, allocatable :: b(:)              ! boundary value
    real, allocatable :: x(:), y(:), z(:)  ! gradient components
    real, allocatable :: q(:)              ! flux of the variable
    real              :: sigma             ! sigma
    real              :: res               ! residual after lin. solver

    ! Boundary cell type (important for scalars, since they
    ! can have different boundary conditions at the walls)
    integer, allocatable :: bnd_cond_type(:)

    ! Parameters for numerical solution of the variable
    character(SL) :: solver          ! solver
    character(SL) :: prec            ! preconditioner
    character(SL) :: o_prec(MAX_STRING_ITEMS)  ! options for preconditioner
    integer       :: adv_scheme      ! advection scheme
    integer       :: grad_method     ! gradient computation method
    real          :: blend           ! blending (1.0 central; 0.0 upwind)
    integer       :: td_scheme       ! time-disretization
    real          :: tol             ! linear solver tolerance
    real          :: urf             ! under-relaxation factor
    integer       :: miter           ! max number of iterations for variable
    integer       :: niter           ! executed number of iterations for var
    logical       :: blend_matrix    ! are you blending matrix with upwind?

    ! Each variable has its own copy of PETSc, how cute is that?
    integer                   :: pet_rank
    type(Petsc_Type), pointer :: Pet
  end type

  contains

#   include "Var_Mod/Create_New_Only.f90"
#   include "Var_Mod/Create_Solution.f90"
#   include "Var_Mod/Destroy_New_Only.f90"
#   include "Var_Mod/Destroy_Solution.f90"
#   include "Var_Mod/Bnd_Cond_Name.f90"
#   include "Var_Mod/Bnd_Cond_Type.f90"

  end module

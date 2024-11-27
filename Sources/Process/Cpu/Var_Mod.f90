#include "../../Shared/Assert.h90"
#include "../../Shared/Unused.h90"

!==============================================================================!
  module Var_Mod
!------------------------------------------------------------------------------!
!>  Module Var_Mod defines Var_Type which is, in turn, used to define all
!>  dependent variables in the T-Flows' solver Process.  Thus, individual
!>  velocity components, temperature, pressure and pressure corrections are all
!>  objects of Var_Type.  The Var_Type holds data fields to store its current
!>  value (n for new) old value (o for old) and older than old (oo) value.
!>  In addition, it also holds fields to store gradients of each dependent
!>  variable as well as many linear solver parameters, from the desired linear
!>  solver, tolerance levels required in linear solvers and even options passed
!>  to PETSc solvers.  Each variable holds its own instance of PETSc solver.
!------------------------------------------------------------------------------!
!   Comments on PETSc integration                                              !
!                                                                              !
!   Var_Mod features unique integration with PETSc object for each variable.   !
!   It maintains an individual PETSc instance (Pet) for optimized performance, !
!   particularly in preconditioning matrix formation and memory management.    !
!   * Advantages:                                                              !
!     - Efficiency: By allocating a dedicated PETSc instance for each          !
!       variable, the module potentially enhances computational efficiency,    !
!       especially in time-intensive operations like preconditioning.          !
!     - Flexibility: The module's design allows for detailed control over each !
!       variable's treatment within the solver, accommodating complex and      !
!       variable-specific requirements.                                        !
!   * Considerations:                                                          !
!     - Memory Usage: The unique PETSc instance per variable approach might    !
!      increase memory consumption, a trade-off for the expected gains in      !
!      computational efficiency.                                               !
!    - Complexity: Managing individual solvers and settings for each variable  !
!      adds complexity to the module, necessitating careful implementation and !
!      documentation.  Furthermore, there is asymmetry with the way in which   !
!      native solvers are treated, namely there is only one instance of them   !
!      in the whole program, as opposed to many PETSc instances.               !
!------------------------------------------------------------------------------!
  use Petsc_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

  !--------------!
  !   Var type   !
  !--------------!
  type Var_Type

    type(Grid_Type),   pointer :: pnt_grid    !! the grid on which it's defined
    type(Matrix_Type), pointer :: pnt_matrix  !! the matrix to store the linear
                                              !! resulting from discretization
    character(VL) :: name       !! variable name, always upper case and
                                !! very short (4, defined in Const_Mod)
    character(VL) :: flux_name  !! variable flux name, always upper case and
                                !! very short (4, defined in Const_Mod)
    real, allocatable :: n(:)   !! new value (current time step)
    real, allocatable :: o(:)   !! old value (previous time step)
    real, allocatable :: oo(:)  !! time step before old (older than old)
    real, allocatable :: b(:)   !! boundary value
    real, allocatable :: x(:)   !! x component of variable's gradient
    real, allocatable :: y(:)   !! y component of variable's gradient
    real, allocatable :: z(:)   !! z component of variable's gradient
    real, allocatable :: q(:)   !! wall flux of the variable
    real              :: sigma  !! sigma, diffusion transport coefficient
    real              :: res    !! residual after linear solver

    ! Boundary cell type (important for scalars, since they
    ! can have different boundary conditions at the walls)
    integer, allocatable :: bnd_cond_type(:)

    ! Parameters for numerical solution of the variable
    character(SL) :: solver        !! linear solver to use for this variable
    character(SL) :: prec          !! preconditioner for this variable
    character(SL) :: o_prec(MAX_STRING_ITEMS)
      !! list of options passed to PETSc preconditioner
    integer       :: adv_scheme    !! advection scheme
    integer       :: grad_method   !! gradient computation method
    real          :: blend         !! upwind blending (1.0 central; 0.0 upwind)
    integer       :: td_scheme     !! time-disretization scheme
    real          :: tol           !! linear solver tolerance
    real          :: urf           !! under-relaxation factor
    integer       :: miter         !! max number of iterations for variable
    integer       :: niter         !! executed number of iterations for var
    logical       :: blend_matrix  !! are you blending matrix with upwind?

    ! Each variable has its own copy of PETSc, how cute is that?
    integer                   :: pet_rank  !! rank of the PETSc object
    type(Petsc_Type), pointer :: Pet       !! PETSc object associated
                                           !! with the variable
  end type

  contains

#   include "Var_Mod/Create_New_Only.f90"
#   include "Var_Mod/Create_Solution.f90"
#   include "Var_Mod/Destroy_New_Only.f90"
#   include "Var_Mod/Destroy_Solution.f90"
#   include "Var_Mod/Bnd_Cond_Name.f90"
#   include "Var_Mod/Bnd_Cond_Type.f90"

  end module

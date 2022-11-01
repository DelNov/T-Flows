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

    type(Grid_Type), pointer :: pnt_grid  ! grid for which it is defined

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
    character(SL) :: prec_opts(MSI)  ! preconditioner options
    integer       :: adv_scheme      ! advection scheme
    integer       :: grad_method     ! gradient computation method
    real          :: blend           ! blending (1.0 central; 0.0 upwind)
    integer       :: td_scheme       ! time-disretization
    real          :: tol             ! linear solver tolerance
    real          :: urf             ! under-relaxation factor
    integer       :: mniter          ! max number of iterations for variable
    integer       :: eniter          ! executed number of iterations for var
  end type

  contains

#   include "Var_Mod/Allocate_New_Only.f90"
#   include "Var_Mod/Allocate_Solution.f90"
#   include "Var_Mod/Bnd_Cond_Name.f90"
#   include "Var_Mod/Bnd_Cond_Type.f90"

  end module

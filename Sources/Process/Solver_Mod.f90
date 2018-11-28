!==============================================================================!
  module Solver_Mod
!------------------------------------------------------------------------------!
!   Module used for native linear solvers.                                     !
!------------------------------------------------------------------------------!
  use Comm_Mod,   only: this_proc
  use Grid_Mod,   only: Grid_Type
  use Matrix_Mod, only: Matrix_Type, Matrix_Mod_Create, Matrix_Mod_Create_Level
  use Vector_Mod, only: Vector_Type, Vector_Mod_Allocate
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !   Solvers type   !
  !------------------!
  type Solver_Type

    type(Grid_Type), pointer :: pnt_grid

    ! Preconditioning "matrix" for single grid methods
    type(Matrix_Type) :: a     ! system matrix for all variables
    type(Matrix_Type) :: d     ! preconditioning "matrix"
    type(Vector_Type) :: b

    ! Hierarchy of linear systems for multigrid methods
    type(Matrix_Type), allocatable :: a_lev(:)  ! system matrix
    type(Matrix_Type), allocatable :: d_lev(:)  ! preconditiong matrix
    type(Vector_Type), allocatable :: x_lev(:)  ! unknown vector
    type(Vector_Type), allocatable :: b_lev(:)  ! right hand side

  end type

  contains

  include 'Solver_Mod/Allocate.f90'
  include 'Solver_Mod/Bicg.f90'
  include 'Solver_Mod/Cg.f90'
  include 'Solver_Mod/Cgs.f90'
  include 'Solver_Mod/Prec_Form.f90'
  include 'Solver_Mod/Prec_Solve.f90'
  include 'Solver_Mod/Residual_Vector.f90'
  include 'Solver_Mod/Normalized_Residual.f90'

  end module 

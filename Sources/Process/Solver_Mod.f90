!==============================================================================!
  module Solver_Mod
!------------------------------------------------------------------------------!
!   Module used for native linear solvers.                                     !
!------------------------------------------------------------------------------!
  use Const_Mod
  use Comm_Mod
  use Grid_Mod,     only: Grid_Type
  use Matrix_Mod
  use Vector_Mod
  use Control_Mod
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
    type(Vector_Type), allocatable :: p_lev(:)  ! one of the CG vectors
    type(Vector_Type), allocatable :: r_lev(:)  ! residual

  end type

  contains

  include 'Solver_Mod/Acm.f90'                  ! additive correction mg
  include 'Solver_Mod/Acm_Coarsen_Matrix.f90'   ! coarsen system matrix
  include 'Solver_Mod/Alias_System.f90'
  include 'Solver_Mod/Bicg.f90'                 ! bicg solver
  include 'Solver_Mod/Cg.f90'                   ! cg solver
  include 'Solver_Mod/Cg_Level.f90'             ! cg smoother for acm
  include 'Solver_Mod/Cgs.f90'                  ! cgs solver
  include 'Solver_Mod/Create.f90'               ! memory all. and creation
  include 'Solver_Mod/Normalized_Root_Mean_Square.f90'
  include 'Solver_Mod/Prec_Form.f90'
  include 'Solver_Mod/Prec_Solve.f90'
  include 'Solver_Mod/Residual_Vector.f90'
  include 'Solver_Mod/Root_Mean_Square.f90'

  end module 

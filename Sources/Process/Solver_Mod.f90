!==============================================================================!
  module Solver_Mod
!------------------------------------------------------------------------------!
!   Module used for native linear solvers.                                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
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

    ! Matrix for all variables except momentum
    type(Matrix_Type) :: a

    ! Matrix for discretized momentum equations
    type(Matrix_Type) :: m

    ! Preconditioning matrix for all variables (used inside solvers only)
    type(Matrix_Type) :: d

    ! Right-hand side for all variables
    ! (used in solvers and during discretization)
    type(Vector_Type) :: b

  end type

  contains

  include 'Solver_Mod/Alias_System.f90'
  include 'Solver_Mod/Bicg.f90'                 ! bicg solver
  include 'Solver_Mod/Cg.f90'                   ! cg solver
  include 'Solver_Mod/Cgs.f90'                  ! cgs solver
  include 'Solver_Mod/Create.f90'               ! memory all. and creation
  include 'Solver_Mod/Normalized_Root_Mean_Square.f90'
  include 'Solver_Mod/Prec_Form.f90'
  include 'Solver_Mod/Prec_Solve.f90'
  include 'Solver_Mod/Residual_Vector.f90'
  include 'Solver_Mod/Root_Mean_Square.f90'

  end module 

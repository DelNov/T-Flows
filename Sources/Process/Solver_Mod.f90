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

    ! Preconditioning "matrix" for single grid methods
    type(Matrix_Type) :: a     ! system matrix for all variables
    type(Matrix_Type) :: d     ! preconditioning "matrix"
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

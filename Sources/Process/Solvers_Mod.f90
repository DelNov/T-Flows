!==============================================================================!
  module Solvers_Mod
!------------------------------------------------------------------------------!
!   Module used for native linear solvers.                                     !
!------------------------------------------------------------------------------!
  use Matrix_Mod
  use Vector_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Preconditioning "matrix" for single grid methods
  type(Matrix_Type) :: d  ! preconditioning "matrix"

  ! Hierarchy of linear systems for multigrid methods
  type(Matrix_Type), allocatable :: a_lev(:)  ! system matrix
  type(Matrix_Type), allocatable :: d_lev(:)  ! preconditiong matrix
  type(Vector_Type), allocatable :: x_lev(:)  ! unknown vector
  type(Vector_Type), allocatable :: b_lev(:)  ! right hand side

  contains

  include 'Solvers_Mod/Bicg.f90'
  include 'Solvers_Mod/Cg.f90'
  include 'Solvers_Mod/Cgs.f90'
  include 'Solvers_Mod/Prec_Form.f90'
  include 'Solvers_Mod/Prec_Solve.f90'
  include 'Solvers_Mod/Residual_Vector.f90'
  include 'Solvers_Mod/Normalized_Residual.f90'

  end module 

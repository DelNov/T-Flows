#include "../Shared/Assert.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Native_Mod
!------------------------------------------------------------------------------!
!>  This module provides functionalities for native (home-grown) linear solvers.
!>  It integrates a few essential components required for solving linear
!>  systems of equations arising from discretization of fluid flow equations.
!>  The module is a bit limited in terms of a number of solvers it offers.
!>  Indeed, it only has Conjugate Gradient (CG) and Bi-Conjugate Gradient
!>  (BiCG) solvers, each with the choice of no preconditioning, diagonal
!>  preconditioning or Incomplete Cholesky preconditioning. It is highly
!>  unlikely that the module will grow further as PETSc libraries, which are
!>  linked with T-Flows, offer a variety of linear solvers and preconditioners.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Matrix_Mod
  use Vector_Mod
  use Control_Mod
  use Work_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------!
  !   Native type   !
  !-----------------!
  !> Encapsulates data structures and procedures for linear solver operations.
  type Native_Type

    type(Grid_Type), pointer :: pnt_grid  !! pointer to the numerical grid

    ! Matrix for all variables except momentum
    type(Matrix_Type) :: A  !! system matrix for all variables except momentum

    ! Matrix for discretized momentum equations
    type(Matrix_Type) :: M  !! system matrix for momentum conservation equations

    ! Right-hand side for all variables
    ! (used in solvers and during discretization)
    type(Vector_Type) :: b  !! right hand side of the linear
                            !! system for all variables
    contains
      procedure, private :: Bicg                 !! biconjugate gradient solver
      procedure, private :: Cg                   !! conjugate gradient solver
      procedure          :: Create_Native        !! creates native solver context
      procedure, private :: Normalized_Root_Mean_Square  !! normalized RMS
                                                         !! computation
      procedure, private :: Prec_Form            !! preconditioner formation
      procedure, private :: Prec_Solve           !! preconditioner solving
      procedure, private :: Residual_Vector      !! computes the residual vector
      procedure, private :: Root_Mean_Square     !! computes root mean square
      procedure          :: Solve_Native         !! primary solver function

  end type

  contains

#   include "Native_Mod/Bicg.f90"
#   include "Native_Mod/Cg.f90"
#   include "Native_Mod/Create_Native.f90"
#   include "Native_Mod/Normalized_Root_Mean_Square.f90"
#   include "Native_Mod/Prec_Form.f90"
#   include "Native_Mod/Prec_Solve.f90"
#   include "Native_Mod/Residual_Vector.f90"
#   include "Native_Mod/Root_Mean_Square.f90"
#   include "Native_Mod/Solve_Native.f90"

  end module

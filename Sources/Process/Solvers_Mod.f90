!==============================================================================!
  module Solvers_Mod
!------------------------------------------------------------------------------!
!   Module used for native linear solvers.                                     !
!------------------------------------------------------------------------------!
  use Matrix_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Preconditioning "matrix" (D)
  type(Matrix_Type) :: D  ! preconditioning "matrix"

  contains

  include 'Solvers_Mod/Bicg.f90'
  include 'Solvers_Mod/Cg.f90'
  include 'Solvers_Mod/Cgs.f90'
  include 'Solvers_Mod/Prec_Form.f90'
  include 'Solvers_Mod/Prec_Solve.f90'
  include 'Solvers_Mod/Residual_Vector.f90'
  include 'Solvers_Mod/Normalized_Residual.f90'

  end module 

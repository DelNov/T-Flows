!==============================================================================!
  subroutine Linear_Solvers(Control, name, verbose)
!------------------------------------------------------------------------------!
!>  Reads which linear solvers to use (native or PETSc) from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  character(SL)       :: name     !! which solvers to use (native or petsc)
  logical, optional   :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Char_Item('LINEAR_SOLVERS',  &
                                'native', name, verbose=verbose)
  call String % To_Upper_Case(name)

  end subroutine

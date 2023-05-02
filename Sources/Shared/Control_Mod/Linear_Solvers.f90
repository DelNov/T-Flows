!==============================================================================!
  subroutine Linear_Solvers(Control, name, verbose)
!------------------------------------------------------------------------------!
!   Reading linear solvers (native or PETSc) from the control file.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  character(SL)       :: name
  logical, optional   :: verbose
!==============================================================================!

  call Control % Read_Char_Item('LINEAR_SOLVERS',  &
                                'native', name, verbose=verbose)
  call String % To_Upper_Case(name)

  end subroutine

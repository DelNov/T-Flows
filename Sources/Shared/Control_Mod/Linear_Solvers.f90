!==============================================================================!
  subroutine Control_Mod_Linear_Solvers(name, verbose)
!------------------------------------------------------------------------------!
!   Reading linear solvers (native or PETSc) from the control file.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL)     :: name  ! name of the pressure-momentum coupling algorithm
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('LINEAR_SOLVERS',  &
                                  'native', name, verbose = .true.)
  call String % To_Upper_Case(name)

  end subroutine

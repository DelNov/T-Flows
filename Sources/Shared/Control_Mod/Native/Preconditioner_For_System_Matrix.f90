!==============================================================================!
  subroutine Preconditioner_For_System_Matrix(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads which preconditioner to use from the file.
!------------------------------------------------------------------------------!
!   These are preconditioners for the T-Flows suite of solvers.  PETSc has     !
!   many more, of course.  For compatibillity with PETSc, these keywords have  !
!   the same values as in PETSc, and all are in lower case.                    !
!                                                                              !
!   It is probably also worth noting that, for native solvers, the only pre-   !
!   conditioner which really makes sense is Incomplete Cholesky.  Jacobi (aka  !
!   diagonal) is slow and no (none) preconditioning makes solvers diverge.     !
!   Having said that, it is sufficient to have one procedure (this) which      !
!   reads preconditioner for all variables, and maybe even that is overkill,   !
!   maybe native solvers should use icc all the time without even asking.      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  character(SL)       :: val      !! preconditioner (diagonal, jacobi, icc, ...)
  logical, optional   :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Char_Item('PRECONDITIONER_FOR_SYSTEM_MATRIX',  &
                                'incomplete_cholesky',               &
                                 val,                                &
                                 verbose)
  call String % To_Lower_Case(val)

  ! Check validity of the input
  if( val .ne. 'none'                .and.  &
      val .ne. 'diagonal'            .and.  &
      val .ne. 'jacobi'              .and.  &
      val .ne. 'incomplete_cholesky' .and.  &
      val .ne. 'icc') then
    call Message % Error(72,                                                &
             'Unknown preconditioner for the system matrix: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',                      &
             file=__FILE__, line=__LINE__, one_proc=.true.)

  ! Use the same names as PETSc
  else
    if(val .eq. 'diagonal')            val = 'jacobi'
    if(val .eq. 'incomplete_cholesky') val = 'icc'
  end if

  end subroutine

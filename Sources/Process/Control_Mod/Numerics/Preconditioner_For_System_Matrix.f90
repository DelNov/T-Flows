!==============================================================================!
  subroutine Control_Mod_Preconditioner_For_System_Matrix(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=80) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('PRECONDITIONER_FOR_SYSTEM_MATRIX',  &
                                  'incomplete_cholesky',               &
                                   val,                                &
                                   verbose)
  call To_Upper_Case(val)

  if( val.ne.'NONE'                .and.  &
      val.ne.'DIAGONAL'            .and.  &
      val.ne.'INCOMPLETE_CHOLESKY' ) then
    print *, '# Unknown preconditioner for the system matrix: ', trim(val)
    print *, '# Exiting!'
    stop 
  end if

  end subroutine

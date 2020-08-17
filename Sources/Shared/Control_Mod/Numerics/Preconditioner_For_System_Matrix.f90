!==============================================================================!
  subroutine Control_Mod_Preconditioner_For_System_Matrix(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL)     :: val
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
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown preconditioner for the system matrix: ',  &
               trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End
  end if

  end subroutine

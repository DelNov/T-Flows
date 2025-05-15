!==============================================================================!
  subroutine Finalize_Work(Work)
!------------------------------------------------------------------------------!
!>  Provides a summary of the usage of working arrays allocated in the Work
!>  object. This subroutine is used for diagnostic and monitoring purposes,
!>  at the end of a program's execution.  It can be used to optimize memory
!>  usage in Work_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type) :: Work
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: T = 50
!-----------------------------------[Locals]-----------------------------------!
  character(DL) :: line
!==============================================================================!

  if(First_Proc()) then
    line( 1:160) = ' '

    write(line(T+1:T+28),'(a)') ' #=========================#'
    print '(a)', trim(line)
    write(line(T+1:T+28),'(a)') ' #     Work_Mod usage      #'
    print '(a)', trim(line)
    write(line(T+1:T+28),'(a)') ' #-------------------------#'
    print '(a)', trim(line)
    write(line(T+1:T+28),'(a,i2,a3)') ' # - Real cell arrays: ',  &
                                      Work % max_r_cell, '  #'
    print '(a)', trim(line)
    write(line(T+1:T+28),'(a,i2,a3)') ' # - Real face arrays: ',  &
                                      Work % max_r_face, '  #'
    print '(a)', trim(line)
    write(line(T+1:T+28),'(a,i2,a3)') ' # - Real node arrays: ',  &
                                      Work % max_r_node, '  #'
    print '(a)', trim(line)
    write(line(T+1:T+28),'(a,i2,a3)') ' # - Int. cell arrays: ',  &
                                      Work % max_i_cell, '  #'
    print '(a)', trim(line)
    write(line(T+1:T+28),'(a,i2,a3)') ' # - Int. face arrays: ',  &
                                      Work % max_i_face, '  #'
    print '(a)', trim(line)
    write(line(T+1:T+28),'(a,i2,a3)') ' # - Int. node arrays: ',  &
                                      Work % max_i_node, '  #'
    print '(a)', trim(line)
    write(line(T+1:T+28),'(a)') ' #-------------------------#'
    print '(a)', trim(line)

  end if

  end subroutine


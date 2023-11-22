!==============================================================================!
  subroutine Finalize_Work(Work)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type) :: Work
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: T = 48
!-----------------------------------[Locals]-----------------------------------!
  character(DL) :: line
!==============================================================================!

  if(First_Proc()) then
    line( 1:160) = ' '

    write(line(T+1:T+33),'(a)') ' #==============================#'
    print '(a)', trim(line)
    write(line(T+1:T+33),'(a)') ' #        Work_Mod usage        #'
    print '(a)', trim(line)
    write(line(T+1:T+33),'(a)') ' #------------------------------#'
    print '(a)', trim(line)
    write(line(T+1:T+33),'(a,i2,a3,i2,a3)') ' # - Real cell arrays: ',  &
                                            Work % max_r_cell, ' / ',   &
                                            Work % req_r_cell, '  #'
    print '(a)', trim(line)
    write(line(T+1:T+33),'(a,i2,a3,i2,a3)') ' # - Real face arrays: ',  &
                                            Work % max_r_face, ' / ',   &
                                            Work % req_r_face, '  #'
    print '(a)', trim(line)
    write(line(T+1:T+33),'(a,i2,a3,i2,a3)') ' # - Real node arrays: ',  &
                                            Work % max_r_node, ' / ',   &
                                            Work % req_r_node, '  #'
    print '(a)', trim(line)
    write(line(T+1:T+33),'(a,i2,a3,i2,a3)') ' # - Int. cell arrays: ',  &
                                            Work % max_i_cell, ' / ',   &
                                            Work % req_i_cell, '  #'
    print '(a)', trim(line)
    write(line(T+1:T+33),'(a,i2,a3,i2,a3)') ' # - Int. face arrays: ',  &
                                            Work % max_i_face, ' / ',   &
                                            Work % req_i_face, '  #'
    print '(a)', trim(line)
    write(line(T+1:T+33),'(a,i2,a3,i2,a3)') ' # - Int. node arrays: ',  &
                                            Work % max_i_node, ' / ',   &
                                            Work % req_i_node, '  #'
    print '(a)', trim(line)
    write(line(T+1:T+33),'(a)') ' #------------------------------#'
    print '(a)', trim(line)

  end if

  end subroutine


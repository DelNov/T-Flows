!==============================================================================!
  subroutine Comm_Mod_Write_Text(fh, text_out, disp)
!------------------------------------------------------------------------------!
!   Write string array for sequential runs.                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer   :: fh            ! file handle
  character :: text_out*(*)  ! text to write out
  integer   :: disp          ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: leng
!==============================================================================!

  leng = len(text_out)

  write(9) text_out

  disp = disp + leng

  end subroutine

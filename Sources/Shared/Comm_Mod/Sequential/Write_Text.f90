!==============================================================================!
  subroutine Write_Text(Comm, fh, text_out, disp)
!------------------------------------------------------------------------------!
!   Write string array for sequential runs                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: fh            ! file handle
  character        :: text_out*(*)  ! text to write out
  integer(DP)      :: disp          ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: leng
!==============================================================================!

  leng = len(text_out)

  write(fh) text_out

  disp = disp + leng

  end subroutine

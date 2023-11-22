!==============================================================================!
  subroutine Write_Text(Comm, fh, text_out, disp)
!------------------------------------------------------------------------------!
!   Write string array for sequential runs                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm
  integer,          intent(in)    :: fh            ! file handle
  character,        intent(in)    :: text_out*(*)  ! text to write out
  integer(DP),      intent(inout) :: disp          ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: leng
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  leng = len(text_out)

  write(fh) text_out

  disp = disp + leng

  end subroutine

!==============================================================================!
  subroutine Read_Text(Comm, fh, text_in, disp)
!------------------------------------------------------------------------------!
!   Read string arrayr for sequential runs.                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm
  integer,          intent(in)    :: fh           ! file handle
  character,        intent(out)   :: text_in*(*)  ! text to read in
  integer(DP),      intent(inout) :: disp         ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: leng
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  leng = len(text_in)

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read string
  read(fh) text_in

  disp = disp + leng

  end subroutine

!==============================================================================!
  subroutine Comm_Mod_Read_Text(fh, text_in, disp)
!------------------------------------------------------------------------------!
!   Read string arrayr for sequential runs.                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer   :: fh           ! file handle
  character :: text_in*(*)  ! text to write out
  integer   :: disp         ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: leng
!==============================================================================!

  leng = len(text_in)

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read string
  read(fh) text_in

  disp = disp + leng

  end subroutine

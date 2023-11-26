!==============================================================================!
  subroutine Write_Text(Comm, fh, text_out, disp)
!------------------------------------------------------------------------------!
!>  Writes a single character string to a file in sequential environment.
!>  This subroutine is used only to write to backup files in Process and is
!>  invoked through Grid's member Comm in Backup_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm          !! communicator from grid
  integer,          intent(in)    :: fh            !! file handle
  character,        intent(in)    :: text_out*(*)  !! text to write out
  integer(DP),      intent(inout) :: disp          !! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: leng
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  leng = len(text_out)

  write(fh) text_out

  disp = disp + leng

  end subroutine

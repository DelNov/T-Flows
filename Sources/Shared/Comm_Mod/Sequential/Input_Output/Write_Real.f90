!==============================================================================!
  subroutine Write_Real(Comm, fh, num, disp)
!------------------------------------------------------------------------------!
!>  Writes a single real number to a file in parallel environment.
!>  This subroutine is used only to write to backup files in Process and is
!>  invoked through Grid's member Comm in Backup_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm  !! communicator from grid
  integer,          intent(in)    :: fh    !! file handle
  real,             intent(in)    :: num   !! number to write
  integer(DP),      intent(inout) :: disp  !! displacement in bytes
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  write(fh) num

  disp = disp + RP

  end subroutine

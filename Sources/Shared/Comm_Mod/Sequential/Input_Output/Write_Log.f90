!==============================================================================!
  subroutine Write_Log(Comm, fh, var, disp)
!------------------------------------------------------------------------------!
!>  Writes a single logical variable to file in sequential environment.
!>  This subroutine is used only to write to backup files in Process and is
!>  invoked through Grid's member Comm in Backup_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm  !! communicator from grid
  integer,          intent(in)    :: fh    !! file handle
  logical,          intent(in)    :: var   !! variable to write
  integer(DP),      intent(inout) :: disp  !! displacement in bytes
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  write(fh) var

  disp = disp + LP

  end subroutine

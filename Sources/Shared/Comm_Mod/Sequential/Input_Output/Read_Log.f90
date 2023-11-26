!==============================================================================!
  subroutine Read_Log(Comm, fh, var, disp)
!------------------------------------------------------------------------------!
!>  Reads a logical value from a file in sequential environment.
!>  This subroutine is used only in Backup_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm  !! communicator from grid
  integer,          intent(in)    :: fh    !! file handle
  logical,          intent(out)   :: var   !! variable to read
  integer(DP),      intent(inout) :: disp  !! displacement in bytes
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read integer
  read(fh) var

  disp = disp + LP

  end subroutine

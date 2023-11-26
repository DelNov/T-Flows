!==============================================================================!
  subroutine Read_Real(Comm, fh, num, disp)
!------------------------------------------------------------------------------!
!>  Reads a real number from a file in sequential environment.
!>  This subroutine is used only in Backup_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm  !! communicator from grid
  integer,          intent(in)    :: fh    !! file handle
  real,             intent(out)   :: num   !! number to read
  integer(DP),      intent(inout) :: disp  !! displacement in bytes
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read integer
  read(fh) num

  disp = disp + RP

  end subroutine

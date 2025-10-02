!==============================================================================!
  subroutine Read_Log_Array(Comm, fh, arr, disp)
!------------------------------------------------------------------------------!
!>  Reads a logical array from a file in sequential environment. This
!>  subroutine is used only to read Swarm data from backup files and is
!>  invoked through the Grid's member Comm in Backup_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm    !! communicator from grid
  integer,          intent(in)    :: fh      !! file handle
  logical,          intent(out)   :: arr(:)  !! array to read
  integer(DP),      intent(inout) :: disp    !! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: length
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  length = size(arr)

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read integer
  read(fh) arr(1:length)

  disp = disp + LP * length

  end subroutine

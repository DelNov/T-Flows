!==============================================================================!
  subroutine Write_Real_Array(Comm, fh, arr, disp)
!------------------------------------------------------------------------------!
!>  Writes a real array to a file in sequential environment. This
!>  subroutine is used only to write Swarm data to backup files and is
!>  invoked through the Grid's member Comm in Backup_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm    !! communicator from grid
  integer,          intent(in)    :: fh      !! file handle
  real,             intent(in)    :: arr(:)  !! array to write
  integer(DP),      intent(inout) :: disp    !! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: length
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  ! Get array's length
  length = size(arr)

  write(fh) arr(1:length)

  disp = disp + RP * length

  end subroutine

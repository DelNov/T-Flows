!==============================================================================!
  subroutine Read_Bnd_Real(Comm, fh, arr, disp)
!------------------------------------------------------------------------------!
!>  Reads a boundary-cell-based (hence, associated with a grid) real array
!>  from a file, in sequential environment. This subroutine is used only to
!>  read from backup files and is invoked through the Grid's member Comm in
!>  Backup_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm    !! communicator from grid
  integer,          intent(in)    :: fh      !! file handle
  real,             intent(out)   :: arr(:)  !! bnd-cell based array to read
  integer(DP),      intent(inout) :: disp    !! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read "distributed" boundary cell array
  do c = 1, Comm % nb_tot
    read(fh) arr(c)
  end do

  disp = disp + Comm % nb_tot * RP

  end subroutine

!==============================================================================!
  subroutine Write_Cell_Real(Comm, fh, array, disp)
!------------------------------------------------------------------------------!
!>  Writes a cell-based (hence, associated with a grid) real array from a
!>  sequential environment to a file.  This subroutine is used only to write
!>  to backup files and is invoked through the Grid's member Comm in Backup_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm      !! communicator from grid
  integer,          intent(in)    :: fh        !! file handle
  real,             intent(in)    :: array(:)  !! bnd-cell-based array to write
  integer(DP),      intent(inout) :: disp      !! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  ! Write "distributed" cell data 
  do c = 1, Comm % nc_tot
    write(fh) array(c)
  end do

  disp = disp + Comm % nc_tot * RP

  end subroutine

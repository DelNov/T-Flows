!==============================================================================!
  subroutine Comm_Mod_Read_Cell_Real(comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Sequential version of reading a "distributed" cell-based array.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Comm_Type) :: comm
  integer         :: fh         ! file handle
  real            :: array(:)
  integer         :: disp       ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read "distributed" cell data 
  do c = 1, comm % nc_tot
    read(fh) array(c)
  end do

  disp = disp + comm % nc_tot * SIZE_REAL

  end subroutine

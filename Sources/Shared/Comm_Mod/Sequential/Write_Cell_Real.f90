!==============================================================================!
  subroutine Comm_Mod_Write_Cell_Real(comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Sequential version of writing a "distributed" cell-based array.            !
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

  ! Write "distributed" cell data 
  do c = 1, comm % nc_tot
    write(9) array(c)
  end do

  disp = disp + comm % nc_tot * SIZE_REAL

  end subroutine

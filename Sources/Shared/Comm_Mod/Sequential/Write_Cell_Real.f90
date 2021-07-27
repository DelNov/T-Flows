!==============================================================================!
  subroutine Write_Cell_Real(Comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Sequential version of writing a "distributed" cell-based array             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: fh         ! file handle
  real             :: array(:)
  integer          :: disp       ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  ! Write "distributed" cell data 
  do c = 1, Comm % nc_tot
    write(fh) array(c)
  end do

  disp = disp + Comm % nc_tot * SIZE_REAL

  end subroutine

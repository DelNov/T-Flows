!==============================================================================!
  subroutine Comm_Mod_Read_Cell_Real(fh, array, disp)
!------------------------------------------------------------------------------!
!   Sequential version of reading a "distributed" cell-based array.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: fh         ! file handle
  real    :: array(:)
  integer :: disp       ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read "distributed" cell data 
  do c = 1, nc_t
    read(fh) array(c)
  end do

  disp = disp + nc_t * SIZE_REAL

  end subroutine

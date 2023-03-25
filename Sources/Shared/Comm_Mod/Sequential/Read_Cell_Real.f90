!==============================================================================!
  subroutine Read_Cell_Real(Comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Sequential version of reading a "distributed" cell-based array.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type),   intent(in)    :: Comm
  integer,            intent(in)    :: fh         ! file handle
  real, dimension(:), intent(out)   :: array(:)   ! array to read
  integer(DP),        intent(inout) :: disp       ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read "distributed" cell data 
  do c = 1, Comm % nc_tot
    read(fh) array(c)
  end do

  disp = disp + Comm % nc_tot * RP

  end subroutine

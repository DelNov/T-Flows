!==============================================================================!
  subroutine Read_Bnd_Real(Comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Sequential version of reading a "distributed" boundary cell-based array.   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: fh         ! file handle
  real             :: array(:)
  integer(DP)      :: disp       ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read "distributed" boundary cell array 
  do c = 1, Comm % nb_tot
    read(fh) array(c)
  end do

  disp = disp + Comm % nb_tot * RP

  end subroutine

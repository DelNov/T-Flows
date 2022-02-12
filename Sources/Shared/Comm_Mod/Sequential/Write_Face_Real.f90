!==============================================================================!
  subroutine Write_Face_Real(Comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Sequential version of writing a "distributed" face-based array             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: fh         ! file handle
  real             :: array(:)
  integer(DP)      :: disp       ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: s
!==============================================================================!

  ! Write "distributed" face data
  do s = 1, Comm % nf_tot
    write(fh) array(s)
  end do

  disp = disp + Comm % nf_tot * SIZE_REAL

  end subroutine

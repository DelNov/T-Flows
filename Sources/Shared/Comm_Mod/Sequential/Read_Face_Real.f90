!==============================================================================!
  subroutine Read_Face_Real(Comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Sequential version of reading a "distributed" face-based array.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: fh         ! file handle
  real             :: array(:)
  integer          :: disp       ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: s
!==============================================================================!

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read "distributed" face data
  do s = 1, Comm % nf_tot
    read(fh) array(s)
  end do

  disp = disp + Comm % nf_tot * SIZE_REAL

  end subroutine

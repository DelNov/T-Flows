!==============================================================================!
  subroutine Comm_Mod_Read_Face_Real(fh, array, disp)
!------------------------------------------------------------------------------!
!   Sequential version of reading a "distributed" face-based array.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: fh         ! file handle
  real    :: array(:)
  integer :: disp       ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: s
!==============================================================================!

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read "distributed" face data 
  do s = 1, nf_t
    read(fh) array(s)
  end do

  disp = disp + nf_t * SIZE_REAL

  end subroutine

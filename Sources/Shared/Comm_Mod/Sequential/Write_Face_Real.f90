!==============================================================================!
  subroutine Comm_Mod_Write_Face_Real(fh, array, disp)
!------------------------------------------------------------------------------!
!   Sequential version of writing a "distributed" face-based array.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: fh         ! file handle
  real    :: array(:)
  integer :: disp       ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: s
!==============================================================================!

  ! Write "distributed" face data 
  do s = 1, nf_t
    write(9) array(s)
  end do

  disp = disp + nf_t * SIZE_REAL

  end subroutine

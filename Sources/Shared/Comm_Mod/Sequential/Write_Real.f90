!==============================================================================!
  subroutine Comm_Mod_Write_Real(fh, num, disp)
!------------------------------------------------------------------------------!
!   Write single real number for sequential runs.                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: fh    ! file handle
  real    :: num   ! number to write out
  integer :: disp  ! displacement in bytes
!==============================================================================!

  write(9) num

  disp = disp + SIZE_REAL

  end subroutine

!==============================================================================!
  subroutine Comm_Mod_Read_Real(fh, num, disp)
!------------------------------------------------------------------------------!
!   Read single real number for sequential runs.                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: fh    ! file handle
  real    :: num   ! number to read
  integer :: disp  ! displacement in bytes
!==============================================================================!

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read integer
  read(fh) num

  disp = disp + SIZE_REAL

  end subroutine

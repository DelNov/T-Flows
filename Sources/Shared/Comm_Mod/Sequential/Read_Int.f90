!==============================================================================!
  subroutine Read_Int(Comm, fh, num, disp)
!------------------------------------------------------------------------------!
!   Read single integer for sequential runs.                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: fh    ! file handle
  integer          :: num   ! number to read
  integer(DP)      :: disp  ! displacement in bytes
!==============================================================================!

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read integer
  read(fh) num

  disp = disp + SIZE_INT

  end subroutine

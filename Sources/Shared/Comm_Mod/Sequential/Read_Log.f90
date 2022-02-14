!==============================================================================!
  subroutine Read_Log(Comm, fh, var, disp)
!------------------------------------------------------------------------------!
!   Read single logical variable for sequential runs                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: fh    ! file handle
  logical          :: var   ! variable to read
  integer(DP)      :: disp  ! displacement in bytes
!==============================================================================!

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read integer
  read(fh) var

  disp = disp + LP

  end subroutine

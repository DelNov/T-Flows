!==============================================================================!
  subroutine Read_Real(Comm, fh, num, disp)
!------------------------------------------------------------------------------!
!   Read single real number for sequential runs.                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: fh    ! file handle
  real             :: num   ! number to read
  integer(DP)      :: disp  ! displacement in bytes
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read integer
  read(fh) num

  disp = disp + RP

  end subroutine

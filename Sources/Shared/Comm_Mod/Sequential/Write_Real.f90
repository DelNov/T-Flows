!==============================================================================!
  subroutine Write_Real(Comm, fh, num, disp)
!------------------------------------------------------------------------------!
!   Write single real number for sequential runs                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: fh    ! file handle
  real             :: num   ! number to write out
  integer          :: disp  ! displacement in bytes
!==============================================================================!

  write(fh) num

  disp = disp + SIZE_REAL

  end subroutine

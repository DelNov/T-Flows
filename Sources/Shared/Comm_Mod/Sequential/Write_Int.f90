!==============================================================================!
  subroutine Write_Int(Comm, fh, num, disp)
!------------------------------------------------------------------------------!
!   Write single integer for sequential runs                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: fh    ! file handle
  integer          :: num   ! number to write out
  integer          :: disp  ! displacement in bytes
!==============================================================================!

  write(fh) num

  disp = disp + SIZE_INT

  end subroutine

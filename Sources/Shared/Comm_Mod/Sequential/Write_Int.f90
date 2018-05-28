!==============================================================================!
  subroutine Comm_Mod_Write_Int(fh, num, disp)
!------------------------------------------------------------------------------!
!   Write single integer for sequential runs.                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: fh    ! file handle
  integer :: num   ! number to write out
  integer :: disp  ! displacement in bytes
!==============================================================================!

  write(9) num

  disp = disp + SIZE_INT

  end subroutine

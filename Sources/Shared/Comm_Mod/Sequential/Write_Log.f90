!==============================================================================!
  subroutine Write_Log(Comm, fh, var, disp)
!------------------------------------------------------------------------------!
!   Write single logical variable for sequential runs                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: fh    ! file handle
  logical          :: var   ! variable to write out
  integer(DP)      :: disp  ! displacement in bytes
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  write(fh) var

  disp = disp + LP

  end subroutine

!==============================================================================!
  subroutine Write_Real(Comm, fh, num, disp)
!------------------------------------------------------------------------------!
!   Write single real number for sequential runs                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm
  integer,          intent(in)    :: fh    ! file handle
  real,             intent(in)    :: num   ! number to write out
  integer(DP),      intent(inout) :: disp  ! displacement in bytes
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  write(fh) num

  disp = disp + RP

  end subroutine

!==============================================================================!
  subroutine Write_Log_Array(Comm, fh, arr, disp)
!------------------------------------------------------------------------------!
!   Write integer array for sequential runs                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type)      :: Comm
  integer               :: fh    ! file handle
  logical, dimension(:) :: arr   ! array to write out
  integer               :: disp  ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: length
!==============================================================================!

  ! Get array's length
  length = size(arr)

  write(fh) arr(1:length)

  disp = disp + SIZE_LOG * length

  end subroutine

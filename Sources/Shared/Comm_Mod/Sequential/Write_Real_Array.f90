!==============================================================================!
  subroutine Write_Real_Array(Comm, fh, arr, disp)
!------------------------------------------------------------------------------!
!   Write real array for sequential runs.                                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type)     :: Comm
  integer              :: fh    ! file handle
  real,   dimension(:) :: arr   ! array to write out
  integer              :: disp  ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: length
!==============================================================================!

  ! Get array's length
  length = size(arr)

  write(9) arr(1:length)

  disp = disp + SIZE_REAL * length

  end subroutine

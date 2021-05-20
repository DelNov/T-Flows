!==============================================================================!
  subroutine Read_Real_Array(Comm, fh, arr, disp)
!------------------------------------------------------------------------------!
!   Read real array for sequential runs.                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type)     :: Comm
  integer              :: fh    ! file handle
  real,   dimension(:) :: arr   ! array to read
  integer              :: disp  ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: length
!==============================================================================!

  length = size(arr)

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read integer
  read(fh) arr(1:length)

  disp = disp + SIZE_REAL * length

  end subroutine

!==============================================================================!
  subroutine Read_Int_Array(Comm, fh, arr, disp)
!------------------------------------------------------------------------------!
!   Read integer array for sequential runs.                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type)      :: Comm
  integer               :: fh    ! file handle
  integer, dimension(:) :: arr   ! array to read
  integer(DP)           :: disp  ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: length
!==============================================================================!

  length = size(arr)

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read integer
  read(fh) arr(1:length)

  disp = disp + IP * length

  end subroutine

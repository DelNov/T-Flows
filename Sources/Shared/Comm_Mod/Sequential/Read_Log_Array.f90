!==============================================================================!
  subroutine Comm_Mod_Read_Log_Array(fh, arr, disp)
!------------------------------------------------------------------------------!
!   Read logical array for sequential runs.                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer               :: fh    ! file handle
  logical, dimension(:) :: arr   ! array to read
  integer               :: disp  ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: length
!==============================================================================!

  length = size(arr)

  ! Position yourself at the right place inside the file
  call fseek(fh, disp, 0)

  ! Read integer
  read(fh) arr(1:length)

  disp = disp + SIZE_LOG * length

  end subroutine

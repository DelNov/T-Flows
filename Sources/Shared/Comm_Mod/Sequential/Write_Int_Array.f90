!==============================================================================!
  subroutine Comm_Mod_Write_Int_Array(fh, arr, disp)
!------------------------------------------------------------------------------!
!   Write integer array for sequential runs.                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer               :: fh    ! file handle
  integer, dimension(:) :: arr   ! array to write out
  integer               :: disp  ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: length
!==============================================================================!

  ! Get array's length
  length = size(arr)

  write(9) arr(1:length)

  disp = disp + SIZE_INT * length

  end subroutine

!==============================================================================!
  subroutine Read_Binary_Real8_Array(File, un, n, reached_end)
!------------------------------------------------------------------------------!
!   Reads an array of double reals from a file in binary format.               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type)  :: File
  integer           :: un  ! unit
  integer           :: n   ! number of items to read
  logical, optional :: reached_end
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  ! If present, assumed the end of file has not been reached
  if(present(reached_end)) then
    reached_end = .false.
  end if

  !--------------------------------------!
  !  Read the number of items you want   !
  !--------------------------------------!
  do i = 1, n
    read(un, end=2) real8_array(i)
  end do

  ! Error trap, if here, you reached the end of file
2 if(present(reached_end)) then
    reached_end = .true.
  end if

  end subroutine

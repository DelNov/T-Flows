!==============================================================================!
  subroutine Read_Binary_Int8_Array(File, un, n, reached_end)
!------------------------------------------------------------------------------!
!>  Reads a specified number of double precision (8-byte) integers from a
!>  binary file into the File % int8_array.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type)  :: File         !! parent class
  integer           :: un           !! file unit
  integer           :: n            !! number of items to read
  logical, optional :: reached_end  !! flag to indicate if the end of the file
                                    !! was reached during the reading process
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(File)
!==============================================================================!

  ! If present, assumed the end of file has not been reached
  if(present(reached_end)) then
    reached_end = .false.
  end if

  !--------------------------------------!
  !  Read the number of items you want   !
  !--------------------------------------!
  do i = 1, n
    read(un, end=2) int8_array(i)
  end do

  ! Error trap, if here, you reached the end of file
2 if(present(reached_end)) then
    reached_end = .true.
  end if

  end subroutine

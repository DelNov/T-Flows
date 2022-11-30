!==============================================================================!
  logical function Is_In_Binary_Format(File, file_name)
!------------------------------------------------------------------------------!
!   Find out if file is written in binary (or ASCII) format.                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type) :: File
  character(len=*) :: file_name
!-----------------------------------[Locals]-----------------------------------!
  integer    :: file_size, file_unit, i
  integer(1) :: byte
!==============================================================================!

  inquire(file=file_name, size=file_size)

  ! Open the file in binary mode
  call File % Open_For_Reading_Binary(file_name, file_unit, verbose=.false.)

  ! Read 32768 bytes and check if any of them is non-ASCII
  do i = 1, min(32768, file_size)  ! don't read beyond the end of file
    read(file_unit) byte
    ! If byte is non-ASCII, the file is not in ASCII format
    if(byte .lt. 0) then   ! > 127 is for one byte signed integer < 0
      close(file_unit)
      Is_In_Binary_Format = .true.
      return
    end if
  end do

  ! All 32768 bytes read were ASCII, assume the file is ASCII too
  close(file_unit)
  Is_In_Binary_Format = .false.

  end function

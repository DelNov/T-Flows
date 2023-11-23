!==============================================================================!
  integer function Line_Length(File, file_unit)
!------------------------------------------------------------------------------!
!>  Calculates the length of a line in a file.  It achieves that by counting
!>  the number of characters in a line up to the newline character.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type) :: File       !! parent class
  integer          :: file_unit  !! unit of the file being analyzed
!-----------------------------------[Locals]-----------------------------------!
  integer    :: length
  integer(1) :: byte
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(File)
!==============================================================================!

  length = 0
  do
    read(file_unit) byte
    if(byte .eq. 10) goto 1  ! end of line
    length = length + 1
  end do
1 continue

  Line_Length = length

  end function

!==============================================================================!
  integer function Line_Length(File, file_unit)
!------------------------------------------------------------------------------!
!   Count how long a line is                                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type) :: File
  integer          :: file_unit
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

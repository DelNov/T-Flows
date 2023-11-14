!==============================================================================!
  pure function First_Upper(String, char_array) result(char_manip)
!------------------------------------------------------------------------------!
!>  Function that converts the first letter of a string to uppercase while the
!>  rest of the string is in lowercase.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(String_Type), intent(in) :: String      !! parent class
  character(len=*),   intent(in) :: char_array  !! input string
  character(:),      allocatable :: char_manip  !! final (manipulated) string
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, value
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(String)
!==============================================================================!

  ! Copy the string from the argument to the local which you will manipulate
  char_manip = trim(char_array)

  ! Then do the necessary manipulations.  First the entire string to lower ...
  do i = 1, len_trim(char_manip)
    value = ichar(char_manip(i:i))
    if(value >= 65 .and. value <= 90) then
      char_manip(i:i) = char(value+32)
    end if
  end do

  ! ... and finally the first character to upper case
  if(len_trim(char_manip) > 0) then
    value = ichar(char_manip(1:1))
    if (value >= 97 .and. value <= 122) then
      char_manip(1:1) = char(value-32)
    end if
  end if

  end function

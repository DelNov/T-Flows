!==============================================================================!
  pure subroutine To_Lower_Case(String, char_array)
!------------------------------------------------------------------------------!
!>  Transforms the entire string to lowercase.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(String_Type), intent(in)    :: String      !! parent class
  character(len=*),   intent(inout) :: char_array  !! string being manipulated
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, value
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(String)
!==============================================================================!

  do i = 1, len_trim(char_array)
    value = ichar(char_array(i:i))
    if(value >= 65 .and. value <= 90) then
      char_array(i:i) = char(value+32)
    end if
  end do

  end subroutine

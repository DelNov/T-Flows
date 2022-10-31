!==============================================================================!
  subroutine To_Upper_Case(String, char_array)
!------------------------------------------------------------------------------!
!   Transforms String to uppercase.                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(String_Type) :: String
  character(len=*)   :: char_array
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, value
!==============================================================================!

  do i = 1, len_trim(char_array)
    value = ichar(char_array(i:i))
    if (value >= 97 .and. value <= 122) then
      char_array(i:i) = char(value-32)
    end if
  end do

  end subroutine

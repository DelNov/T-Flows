!==============================================================================!
  integer function Key_Ind(key, key_to_values, n_keys)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer      :: n_keys
  character(*) :: key
  character(*) :: key_to_values(n_keys)
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  Key_Ind = 0

  do i = 1, n_keys
    if(key == key_to_values(i)) then
      Key_Ind = i
      return
    end if
  end do

  end function

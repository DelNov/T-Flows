!==============================================================================!
  integer function Key_Ind(key, key_to_values, n_keys)
!------------------------------------------------------------------------------!
!>  Simple lookup function that searches for a specified key within an array
!>  of keys and returns the index of the first occurrence of this key. If the
!>  key is not found, the function returns zero.  It is primarily used when
!>  reading control file, searching for boundary and initial conditions.
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
    if(key .eq. key_to_values(i)) then
      Key_Ind = i
      return
    end if
  end do

  end function

!==============================================================================!
  subroutine Convert_Problem_Name_To_Integer(problem_name, int_p_name)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(*) :: problem_name
  integer      :: int_p_name
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, l, ii
!==============================================================================!

  l = len_trim(problem_name)
  int_p_name = 0

  do ii = 1, l
    i = ichar(problem_name(ii:ii))
    if (mod(ii,2) .ne. 0) then
      int_p_name = int_p_name + i * 1000
    else
      int_p_name = int_p_name + i
    end if
  end do

  end subroutine


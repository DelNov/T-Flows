!==============================================================================!
  subroutine Work_Out_I_Ranges(Mem, i, i_range, i_inc,  &
                                    i_lower, i_upper,   &
                                    file, line)
!------------------------------------------------------------------------------!
!>  Work out the ranges for rows (first index, index i)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Memory_Type), intent(in)  :: Mem         !! parent class
  integer, optional,  intent(in)  :: i           !! matrix index
  integer, optional,  intent(in)  :: i_range(2)  !! range in i
  integer, optional,  intent(in)  :: i_inc       !! size increment
  integer,            intent(out) :: i_lower     !! lower bound of i
  integer,            intent(out) :: i_upper     !! upper bound of i
  character(*),       intent(in)  :: file
  integer,            intent(in)  :: line
!==============================================================================!

  !-------------------------------------------!
  !   Handle validity of optional arguments   !
  !     and work out the i ranges (rows)      !
  !-------------------------------------------!

  ! Arguments i and i_range
  if(.not. present(i_range)) then
    if(present(i)) then
      i_lower = i
      i_upper = i
    else
      ! Neither i, nor the i_range is present
      call Message % Error(48,                                               &
                      'Wrong invocation, either i or i_range should be  '//  &
                      'specified!  This error is critical.  Exiting!',       &
                      file=file, line=line)
    end if
  else
    if(.not. present(i)) then
      i_lower = i_range(1)
      i_upper = i_range(2)
    else
      ! Both i and i_range are present
      call Message % Error(60,                                               &
                      'Wrong invocation, either i or i_range should be  '//  &
                      'specified, but not both of them!  '               //  &
                      'This error is critical.  Exiting!',                   &
                      file=file, line=line)
    end if
  end if

  Assert(i_lower .le. i_upper)

  end subroutine

!==============================================================================!
  subroutine Work_Out_J_Ranges(Mem, j, j_range, j_inc,  &
                                    j_lower, j_upper,   &
                                    file, line)
!------------------------------------------------------------------------------!
!>  Work out the ranges for columns (second index, index j)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Memory_Type), intent(in)  :: Mem         !! parent class
  integer, optional,  intent(in)  :: j           !! matrix index
  integer, optional,  intent(in)  :: j_range(2)  !! matrix range in j
  integer, optional,  intent(in)  :: j_inc       !! size increment
  integer,            intent(out) :: j_lower     !! lower bound of j
  integer,            intent(out) :: j_upper     !! upper bound of j
  character(*),       intent(in)  :: file
  integer,            intent(in)  :: line
!==============================================================================!

  !-------------------------------------------!
  !   Handle validity of optional arguments   !
  !    and work out the j ranges (columns)    !
  !-------------------------------------------!

  if(.not. present(j_range)) then
    if(present(j)) then
      j_lower = j
      j_upper = j
    else
      call Message % Error(48,                                               &
                      'Wrong invocation, either j or j_range shuold be  '//  &
                      'specified!  This error is critical.  Exiting!',       &
                      file=file, line=line)
    end if
  else
    if(.not. present(j)) then
      j_lower = j_range(1)
      j_upper = j_range(2)
    else
      ! Both j and j_range are present
      call Message % Error(60,                                               &
                      'Wrong invocation, either j or j_range should be  '//  &
                      'specified, but not both of them!  '               //  &
                      'This error is critical.  Exiting!',                   &
                      file=file, line=line)
    end if
  end if

  Assert(j_lower .le. j_upper)

  end subroutine

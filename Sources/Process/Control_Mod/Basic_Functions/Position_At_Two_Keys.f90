!==============================================================================!
  subroutine Control_Mod_Position_At_Two_Keys(keyword_1, keyword_2,  &
                                              found, verbose)
!------------------------------------------------------------------------------!
!   Position yourself within the file at the line specified with two keys.     !
!   It is intended to be used to find the boundary condition specifications.   !
!------------------------------------------------------------------------------!
  use Const_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*)  :: keyword_1
  character(len=*)  :: keyword_2
  integer           :: found
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  logical :: reached_end
!==============================================================================!

  rewind(CONTROL_FILE_UNIT)

  ! Browse through command file to find two keywords in one file
  found = NO
  do
    call Tokenizer_Mod_Read_Line(CONTROL_FILE_UNIT, reached_end)
    if(reached_end) goto 1

    call To_Upper_Case(line % tokens(2))

    if(line % tokens(1) .eq. trim(keyword_1) .and.  &
       line % tokens(2) .eq. trim(keyword_2)) then
      found = YES
      return 
    end if
  end do

1 if(found .eq. NO) then
    if(present(verbose)) then
      if(verbose) then
        print '(5a)', ' # Could not find the line with keywords: ',  &
                      keyword_1, ', ', trim(keyword_2), '!'
      end if
    end if
  end if 

  end subroutine

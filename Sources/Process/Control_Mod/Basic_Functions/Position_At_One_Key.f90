!==============================================================================!
  subroutine Control_Mod_Position_At_One_Key(keyword, found, verbose)
!------------------------------------------------------------------------------!
!   Position yourself within the file at the line specified with one key.      !
!   It is intended to be used to find the initial condition specifications.    !
!------------------------------------------------------------------------------!
  use Const_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*)  :: keyword
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

    if(line % tokens(1) .eq. trim(keyword)) then
      found = YES
      return 
    end if
  end do

1 if(found .eq. NO) then
    if(present(verbose)) then
      if(verbose) then
        print '(3a)', ' # Could not find the line with keyword: ',  &
                      keyword, '!'
      end if
    end if
  end if 

  end subroutine

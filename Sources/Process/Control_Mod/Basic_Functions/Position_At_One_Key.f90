!==============================================================================!
  subroutine Control_Mod_Position_At_One_Key(keyword, found, verbose)
!------------------------------------------------------------------------------!
!   Position yourself within the file at the line specified with one key.      !
!   It is intended to be used to find the initial condition specifications.    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*)  :: keyword
  logical           :: found
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  logical :: reached_end
!==============================================================================!

  rewind(CONTROL_FILE_UNIT)

  !-----------------------------------------------------!
  !   Browse through command file to find one keyword   !
  !-----------------------------------------------------!
  found = .false.
  do
    call Tokenizer_Mod_Read_Line(CONTROL_FILE_UNIT, reached_end)
    if(reached_end) goto 1

    ! Found the correct keyword
    if(line % tokens(1) .eq. trim(keyword)) then
      found = .true.
      return

    ! Keyword not found, try to see if there is similar, maybe it was a typo
    else
      call Control_Mod_Similar_Warning(trim(keyword),           &
                                       trim(line % tokens(1)))
    end if

  end do

  !--------------------------------------------!
  !   Keyword was not found; issue a warning   !
  !--------------------------------------------!
1 if(.not. found) then
    if(present(verbose)) then
      if(verbose .and. this_proc < 2) then
        print '(3a)', ' # WARING!  Could not find the line with keyword: ',  &
                      keyword, '!'
      end if
    end if
  end if 

  end subroutine

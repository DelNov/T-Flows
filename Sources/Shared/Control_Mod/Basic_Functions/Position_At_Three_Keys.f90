!==============================================================================!
  subroutine Control_Mod_Position_At_Three_Keys(keyword_1,  &
                                                keyword_2,  &
                                                keyword_3,  &
                                                found, verbose)
!------------------------------------------------------------------------------!
!   Position yourself within the file at the line specified with three keys.   !
!   It is intended to be used to find the interface condition specifications.  !
!                                                                              !
!   This function is case sensitive!                                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*)  :: keyword_1
  character(len=*)  :: keyword_2
  character(len=*)  :: keyword_3
  logical           :: found
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  logical :: reached_end
!==============================================================================!

  rewind(control_file_unit)

  !------------------------------------------------------------------!
  !   Browse through command file to find two keywords in one file   !
  !------------------------------------------------------------------!
  found = .false.
  do
    call File_Mod_Read_Line(control_file_unit, reached_end)
    if(reached_end) goto 1

    ! First keyword is "INTERFACE_CONDITION", ...
    ! ... second and third are two problem names
    if(line % tokens(1) .eq. trim(keyword_1) .and.  &
       line % tokens(2) .eq. trim(keyword_2) .and.  &
       line % tokens(3) .eq. trim(keyword_3)) then
      found = .true.
      return

    end if
  end do

  !--------------------------------------------!
  !   Keyword was not found; issue a warning   !
  !--------------------------------------------!
1 if(.not. found) then
    if(present(verbose)) then
      if(verbose .and. this_proc < 2) then
        print '(5a)', ' # NOTE! Could not find the line with keywords: ',  &
                      keyword_1, ', ', trim(keyword_2), ', ', trim(keyword_3), '!'
      end if
    end if
  end if

  end subroutine

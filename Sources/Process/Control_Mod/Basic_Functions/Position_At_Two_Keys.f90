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
  logical           :: found
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  logical :: reached_end
!==============================================================================!

  rewind(CONTROL_FILE_UNIT)

  !------------------------------------------------------------------!
  !   Browse through command file to find two keywords in one file   !
  !------------------------------------------------------------------!
  found = .false.
  do
    call Tokenizer_Mod_Read_Line(CONTROL_FILE_UNIT, reached_end)
    if(reached_end) goto 1

    call To_Upper_Case(line % tokens(2))

    ! First keyword is "BOUNDARY_CONDITION", ...
    ! ... second is boundary condition name
    if(line % tokens(1) .eq. trim(keyword_1) .and.  &
       line % tokens(2) .eq. trim(keyword_2)) then
      found = .true.
      return

    ! Keywords not found, try to see if there is similar, maybe there was a typo
    else
      call Control_Mod_Similar_Warning(trim(keyword_1),         &
                                       trim(line % tokens(1)))
      call Control_Mod_Similar_Warning(trim(keyword_2),         &
                                       trim(line % tokens(2)),  &
                                       key_type='boundary condition')
    end if
  end do

  !--------------------------------------------!
  !   Keyword was not found; issue a warning   !
  !--------------------------------------------!
1 if(.not. found) then
    if(present(verbose)) then
      if(verbose) then
        print '(5a)', ' # NOTE! Could not find the line with keywords: ',  &
                      keyword_1, ', ', trim(keyword_2), '!'
      end if
    end if
  end if 

  end subroutine

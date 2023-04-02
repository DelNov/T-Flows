!==============================================================================!
  subroutine Position_At_Three_Keys(Control,    &
                                    keyword_1,  &
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
  class(Control_Type) :: Control
  character(len=*)    :: keyword_1
  character(len=*)    :: keyword_2
  character(len=*)    :: keyword_3
  logical             :: found
  logical,   optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  logical :: reached_end
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Control)
!==============================================================================!

  rewind(control_file_unit)

  !------------------------------------------------------------------!
  !   Browse through command file to find two keywords in one file   !
  !------------------------------------------------------------------!
  found = .false.
  do
    call File % Read_Line(control_file_unit, reached_end)
    if(reached_end) goto 1

    ! First keyword is "INTERFACE_CONDITION", ...
    ! ... second and third are two problem names
    if(Line % tokens(1) .eq. trim(keyword_1) .and.  &
       Line % tokens(2) .eq. trim(keyword_2) .and.  &
       Line % tokens(3) .eq. trim(keyword_3)) then
      found = .true.
      return
    end if

  end do

  !--------------------------------------------!
  !   Keyword was not found; issue a warning   !
  !--------------------------------------------!
1 if(.not. found) then
    if(present(verbose)) then
      if(verbose .and. First_Proc()) then
        print '(5a)', ' # NOTE! Could not find the line with keywords: ',  &
                      keyword_1, ', ', trim(keyword_2), ', ', trim(keyword_3), '!'
      end if
    end if
  end if

  end subroutine

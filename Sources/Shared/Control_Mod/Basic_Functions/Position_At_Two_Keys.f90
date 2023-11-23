!==============================================================================!
  subroutine Position_At_Two_Keys(Control,    &
                                  keyword_1,  &
                                  keyword_2,  &
                                  found, verbose)
!------------------------------------------------------------------------------!
!>  Position yourself within the file at the line specified with two keys.
!>  It is intended to be used to find the boundary condition specifications.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)              :: Control    !! parent class
  character(len=*),    intent(in)  :: keyword_1  !! keyword 1
  character(len=*),    intent(in)  :: keyword_2  !! keyword 2
  logical,             intent(out) :: found      !! true if two keys are found
  logical,   optional, intent(in)  :: verbose    !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  logical :: reached_end
!==============================================================================!

  rewind(Control % file_unit)

  !------------------------------------------------------------------!
  !   Browse through command file to find two keywords in one file   !
  !------------------------------------------------------------------!
  found = .false.
  do
    call File % Read_Line(Control % file_unit, reached_end)
    if(reached_end) goto 1

    call String % To_Upper_Case(Line % tokens(2))

    ! First keyword is "BOUNDARY_CONDITION", ...
    ! ... second is boundary condition name
    if(Line % tokens(1) .eq. trim(keyword_1) .and.  &
       Line % tokens(2) .eq. trim(keyword_2)) then
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
                      keyword_1, ', ', trim(keyword_2), '!'
      end if
    end if
  end if

  end subroutine

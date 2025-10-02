!==============================================================================!
  subroutine Position_At_One_Key(Control, keyword, found, verbose)
!------------------------------------------------------------------------------!
!>  Position yourself within the contrl file at the line specified with one key.
!>  It is intended to be used to find the initial condition specifications.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)              :: Control  !! parent class
  character(len=*),    intent(in)  :: keyword  !! keyword it searches
  logical,             intent(out) :: found    !! true if keyword is found
  logical, optional,   intent(in)  :: verbose  !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  logical :: reached_end
!==============================================================================!

  rewind(Control % file_unit)

  !-----------------------------------------------------!
  !   Browse through command file to find one keyword   !
  !-----------------------------------------------------!
  found = .false.
  do
    call File % Read_Line(Control % file_unit, reached_end)
    if(reached_end) goto 1

    ! Found the correct keyword
    if(Line % tokens(1) .eq. trim(keyword)) then
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
        print '(3a)', ' # WARING!  Could not find the line with keyword: ',  &
                      keyword, '!'
      end if
    end if
  end if 

  end subroutine

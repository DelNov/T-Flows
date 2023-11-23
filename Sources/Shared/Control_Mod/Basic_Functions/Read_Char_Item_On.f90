!==============================================================================!
  subroutine Read_Char_Item_On(Control, keyword, def, val, verbose)
!------------------------------------------------------------------------------!
!>  Working horse function to read string values (argument "val") behind a
!>  keyword (argument "keyword") in the control file, staring from the current
!>  position. If not found, a default value specified in argument "def" is used.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)              :: Control  !! parent class
  character(len=*),    intent(in)  :: keyword  !! keyword it searches
  character(len=*),    intent(in)  :: def      !! default value
  character(SL),       intent(out) :: val      !! spefified value, if found
  logical,   optional, intent(in)  :: verbose  !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  logical :: reached_end
!==============================================================================!

  ! Set default value
  val = def

  !---------------------------------------------------------!
  !   Read one line from command file to find the keyword   !
  !---------------------------------------------------------!
  call File % Read_Line(Control % file_unit, reached_end)
  if(reached_end) goto 1

  ! Found the correct keyword
  if(Line % tokens(1) .eq. trim(keyword)) then
    read(Line % tokens(2), *) val
    return
  end if

  !--------------------------------------------!
  !   Keyword was not found; issue a warning   !
  !--------------------------------------------!
1 if(present(verbose)) then
    if(verbose .and. First_Proc()) then
      print '(4a)', ' # NOTE! Could not find the keyword: ',  &
                      trim(keyword), '. Using the default: ', trim(def)
    end if
  end if

  end subroutine

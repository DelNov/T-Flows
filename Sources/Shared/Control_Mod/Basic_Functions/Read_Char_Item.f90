!==============================================================================!
  subroutine Control_Mod_Read_Char_Item(keyword, def, val, verbose)
!------------------------------------------------------------------------------!
!   Working horse function to read integer value (argument "val") behind a     !
!   keyword (argument "keyword") in control file.  If not found, a default     !
!   value specified in argument "def" is used.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*)              :: keyword
  character(len=*), intent(in)  :: def      ! default value
  character(SL),    intent(out) :: val      ! spefified value, if found
  logical,          optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  logical :: reached_end
!==============================================================================!

  rewind(control_file_unit)

  ! Set default value
  val = def

  !--------------------------------------------------------------!
  !   Browse through command file to see if user specificed it   !
  !--------------------------------------------------------------!
  do
    call File % Read_Line(control_file_unit, reached_end)
    if(reached_end) goto 1

    ! Found the correct keyword
    if(Line % tokens(1) .eq. trim(keyword)) then
      read(Line % tokens(2), *) val
      return

    ! Keyword not found, try to see if there is similar, maybe it was a typo
    else
      call Control_Mod_Similar_Warning( keyword, trim(Line % tokens(1)) )
    end if

  end do

  !--------------------------------------------!
  !   Keyword was not found; issue a warning   !
  !--------------------------------------------!
1 if(present(verbose)) then
    if(verbose .and. this_proc < 2) then
      print '(4a)', ' # NOTE! Could not find the keyword: ',  &
                      trim(keyword), '. Using the default: ', trim(def)
    end if
  end if

  end subroutine

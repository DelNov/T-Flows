!==============================================================================!
  subroutine Control_Mod_Read_Char_Item_On(keyword, def, val, verbose)
!------------------------------------------------------------------------------!
!   Working horse function to read integer value (argument "val") behind a     !
!   keyword (argument "keyword") in control file.  If not found, a default     !
!   vaue specified in argument "def" is used.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*)  :: keyword
  character(len=*)  :: def      ! default value
  character(SL)     :: val      ! spefified value, if found
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  logical :: reached_end
!==============================================================================!

  ! Set default value
  val = def

  !---------------------------------------------------------!
  !   Read one line from command file to find the keyword   !
  !---------------------------------------------------------!
  call File_Mod_Read_Line(control_file_unit, reached_end)
  if(reached_end) goto 1

  ! Found the correct keyword
  if(line % tokens(1) .eq. trim(keyword)) then
    read(line % tokens(2), *) val
    return

  ! Keyword not found, try to see if there is similar, maybe it was a typo
  else
    call Control_Mod_Similar_Warning( keyword, trim(line % tokens(1)) )
  end if

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

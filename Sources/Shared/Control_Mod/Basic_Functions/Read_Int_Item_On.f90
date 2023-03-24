!==============================================================================!
  subroutine Control_Mod_Read_Int_Item_On(keyword, def, val, verbose)
!------------------------------------------------------------------------------!
!   Working horse function to read integer value (argument "val") behind a     !
!   keyword (argument "keyword") in control file.  If not found, a default     !
!   vaue specified in argument "def" is used.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*)     :: keyword
  integer, intent(in)  :: def      ! default value
  integer, intent(out) :: val      ! spefified value, if found
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  logical :: reached_end
!==============================================================================!

  ! Set default value
  val = def

  !-----------------------------------------------------!
  !   Browse through command file to find the keyword   !
  !-----------------------------------------------------!
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

  !--------------------------------------------!
  !   Keyword was not found; issue a warning   !
  !--------------------------------------------!
1 if(present(verbose)) then
    if(verbose .and. First_Proc()) then
      if(def .eq. HUGE_INT) then
        print '(3a,i9)', ' # NOTE! Could not find the keyword: ',  &
                          trim(keyword), '. Using the default HUGE_INT value'
      else
        print '(3a,i9)', ' # NOTE! Could not find the keyword: ',  &
                          trim(keyword), '. Using the default: ', def
      end if
    end if
  end if

  end subroutine

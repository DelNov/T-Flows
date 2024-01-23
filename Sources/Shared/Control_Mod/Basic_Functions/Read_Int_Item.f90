!==============================================================================!
  subroutine Read_Int_Item(Control, keyword, def, val, verbose)
!------------------------------------------------------------------------------!
!   Working horse function to read integer value (argument "val") behind a     !
!   keyword (argument "keyword") in control file.  If not found, a default     !
!   value specified in argument "def" is used.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)              :: Control
  character(len=*),    intent(in)  :: keyword
  integer,             intent(in)  :: def      ! default value
  integer,             intent(out) :: val      ! spefified value, if found
  logical,   optional, intent(in)  :: verbose
!-----------------------------------[Locals]-----------------------------------!
  logical :: reached_end
!==============================================================================!

  rewind(Control % file_unit)

  ! Set default value
  val = def

  !-----------------------------------------------------!
  !   Browse through command file to find the keyword   !
  !-----------------------------------------------------!
  do
    call File % Read_Line(Control % file_unit, reached_end)
    if(reached_end) goto 1

    ! Found the correct keyword
    if(Line % tokens(1) .eq. trim(keyword)) then
      read(Line % tokens(2), *) val
      return
    end if

  end do

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

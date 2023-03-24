!==============================================================================!
  subroutine Control_Mod_Read_Real_Array_On(keyword, values, n, verbose)
!------------------------------------------------------------------------------!
!   Working horse function to read integer value (argument "val") behind a     !
!   keyword (argument "keyword") in control file.  If not found, a default     !
!   vaue specified in argument "def" is used.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*)  :: keyword
  real              :: values(128)   ! spefified value, if found
  integer           :: n             ! number of items
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
  logical :: reached_end
!==============================================================================!

  ! Set default values
  values = 0.0

  !---------------------------------------------------------!
  !   Read one line from command file to find the keyword   !
  !---------------------------------------------------------!
  call File % Read_Line(control_file_unit, reached_end)
  if(reached_end) goto 1

  ! Found the correct keyword
  if(Line % tokens(1) .eq. trim(keyword)) then

    do i=2, Line % n_tokens
      read(Line % tokens(i), *) values(i-1)
    end do
    n = Line % n_tokens - 1
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
      print '(2a)', ' # NOTE! Could not find the keyword: ',  &
                      trim(keyword)
    end if
  end if
  n = 0

  end subroutine

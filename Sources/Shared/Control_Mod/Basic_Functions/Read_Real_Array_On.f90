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
  call File_Mod_Read_Line(control_file_unit, reached_end)
  if(reached_end) goto 1

  ! Found the correct keyword
  if(line % tokens(1) .eq. trim(keyword)) then

    do i=2, line % n_tokens
      read(line % tokens(i), *) values(i-1)
    end do
    n = line % n_tokens - 1
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
      print '(2a)', ' # NOTE! Could not find the keyword: ',  &
                      trim(keyword)
    end if
  end if
  n = 0

  end subroutine
